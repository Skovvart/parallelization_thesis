module skovvart.actulus.cubase.Device

open Microsoft.FSharp.Quotations
open Alea.CUDA
open Alea.CUDA.Utilities
open skovvart.actulus.cubase.Common

//Reflected definitions for use in quotations. 
[<ReflectedDefinition>]
let arrayCopy (src:deviceptr<'T>) (destination:deviceptr<'T>) length =
    for i = 0 to length - 1 do
        destination.[i] <- src.[i]

//inline for genericism
[<ReflectedDefinition>]
let inline sax a (x:deviceptr<'T>) (res:deviceptr<'T>) length =
    for i = 0 to length - 1 do
        res.[i] <- a * x.[i]

[<ReflectedDefinition>]
let inline saxpy a (x:deviceptr<'T>) (y:deviceptr<'T>) (res:deviceptr<'T>) length =
    for i = 0 to length - 1 do
        res.[i] <- a * x.[i] + y.[i]

let inline RK4_n (plan:Plan<_>) = cuda {
    let states = plan.states
    let! kernel =
        <@ fun (a:int) b steps (Va:deviceptr<'T>) (d_params:deviceptr<'T>) (result:deviceptr<'T>) ->
            let offset = (blockIdx.x * blockDim.x + threadIdx.x) * states * (a + 1)

//not enough shared memory per block to hold all params (no point either, I guess), block-specific params seems slower
//could experiment with just having the thread-specific params copied to shared memory, but would it be worth it for 4 dV/bj_ii calls?
//answer - seems not, the access-time doesn't make up for the copying from global memory, at least not with 2 params

//Block shared memory
//            let blockParamSize = plan.paramCount*1024
//            let paramOffset = blockIdx.x*blockDim.x*plan.paramCount
//            let sharedParams = __shared__.Array<'T> blockParamSize |> __array_to_ptr
//            if(threadIdx.x = 0) then
//                arrayCopy (d_params.Ptr paramOffset) sharedParams blockParamSize
//            __syncthreads() 
//            let localParams = sharedParams.Ptr (threadIdx.x*plan.paramCount)
//end block shared memory

//Thread local
//            let localParams = __local__.Array plan.paramCount |> __array_to_ptr
//            arrayCopy (d_params.Ptr ((threadIdx.x + blockIdx.x*blockDim.x)*plan.paramCount)) localParams plan.paramCount
//end thread local

//Thread local
//            let localParams = __shared__.Array plan.paramCount |> __array_to_ptr
//            arrayCopy (d_params.Ptr ((threadIdx.x + blockIdx.x*blockDim.x)*plan.paramCount)) localParams plan.paramCount
//end thread local

//global
            let localParams = d_params.Ptr ((threadIdx.x + blockIdx.x*blockDim.x)*plan.paramCount)
//end global    
            let dV = %plan.dV
            let bj_ii = %plan.bj_ii
            let intConv = %intConv
            let floatConv = %floatConv

            let h   = (floatConv -1.0) / intConv steps
            let k1  = __local__.Array<'T>(states) |> __array_to_ptr
            let k2  = __local__.Array<'T>(states) |> __array_to_ptr
            let k3  = __local__.Array<'T>(states) |> __array_to_ptr
            let k4  = __local__.Array<'T>(states) |> __array_to_ptr
            let tmp = __local__.Array<'T>(states) |> __array_to_ptr
            let v   = __local__.Array<'T>(states) |> __array_to_ptr
            
            arrayCopy Va (result.Ptr (offset + a * states)) states
            let mutable y = a
            while y > b do
                bj_ii (intConv y) localParams v
                saxpy (intConv 1) v (result.Ptr (offset + y * states)) v states
                for s = 0 to steps - 1 do // Integrate backwards over [y, y-1]
                    let t = intConv y - intConv s / intConv steps
                    // Hack: Fake limit from left
                    dV (if s = 0 then t - floatConv 1e-5 else t) v localParams k1
                    sax h k1 k1 states
                    saxpy (floatConv 0.5) k1 v tmp states
                    dV (t + h / intConv 2) tmp localParams k2
                    sax h k2 k2 states
                    saxpy (floatConv 0.5) k2 v tmp states
                    dV (t + h / intConv 2) tmp localParams k3
                    sax h k3 k3 states
                    saxpy (intConv 1) k3 v tmp states
                    // Hack: Fake limit from right
                    dV (if s = steps - 1 then t + h + floatConv 1e-5 else t + h) tmp localParams k4
                    sax h k4 k4 states
                    saxpy (intConv 1 / intConv 6) k4 v tmp states
                    saxpy (intConv 1 / intConv 3) k3 tmp tmp states
                    saxpy (intConv 1 / intConv 3) k2 tmp tmp states
                    saxpy (intConv 1 / intConv 6) k1 tmp v states
                arrayCopy v (result.Ptr (offset + y*states - states)) states
                y <- y - 1
        @> |> Compiler.DefineKernel 

    return Entry(fun program ->
        let worker = program.Worker
        let kernel = program.Apply kernel
//        printfn "Kernel{BinaryVersion:%d; NumRegs:%d; PTXVersion:%d; SharedSizeBytes:%d; ConstSizeBytes:%d; LocalSizeBytes:%d;}" 
//            kernel.BinaryVersion kernel.NumRegs kernel.PTXVersion kernel.SharedSizeBytes kernel.ConstSizeBytes kernel.LocalSizeBytes 
        
        fun a b steps blockSize gridSize ->
            let n = (a + 1) * states * blockSize * gridSize
            use Va = worker.Malloc<'T>(Array.zeroCreate states)
            use planParams = worker.Malloc<'T>(plan.planParams)
            use result = worker.Malloc<'T>(Array.zeroCreate n)

            let lp = LaunchParam (gridSize, blockSize)
            let msec, _ = time false "" (fun () -> kernel.Launch lp a b steps Va.Ptr planParams.Ptr result.Ptr)
            let result = result.Gather()

            result, msec
        )
}