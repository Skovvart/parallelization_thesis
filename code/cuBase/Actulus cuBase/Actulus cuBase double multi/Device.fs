module skovvart.actulus.cubase.Device

open Microsoft.FSharp.Quotations
open Alea.CUDA
open Alea.CUDA.Utilities

[<ReflectedDefinition>]
let arrayCopy (src:deviceptr<float>) (destination:deviceptr<float>) length =
    for i = 0 to length - 1 do
        destination.[i] <- src.[i]

[<ReflectedDefinition>]
let sax a (x:deviceptr<float>) (res:deviceptr<float>) length =
    for i = 0 to length - 1 do
        res.[i] <- a * x.[i]

[<ReflectedDefinition>]
let saxpy a (x:deviceptr<float>) (y:deviceptr<float>) (res:deviceptr<float>) length =
    for i = 0 to length - 1 do
        res.[i] <- a * x.[i] + y.[i]

let RK4_n dV bj_ii states = cuda {
    let! kernel =
        <@ fun (a:int) b steps (Va:deviceptr<float>) (result:deviceptr<float>) ->
            let offset = (blockIdx.x * blockDim.x + threadIdx.x) * states * (a + 1)
            let dV = %dV
            let bj_ii = %bj_ii

            let h   = -1.0 / float steps
            let k1  = __local__.Array<float>(states) |> __array_to_ptr
            let k2  = __local__.Array<float>(states) |> __array_to_ptr
            let k3  = __local__.Array<float>(states) |> __array_to_ptr
            let k4  = __local__.Array<float>(states) |> __array_to_ptr
            let tmp = __local__.Array<float>(states) |> __array_to_ptr
            let v   = __local__.Array<float>(states) |> __array_to_ptr
            
            arrayCopy Va (result.Ptr (offset + a * states)) states
            let mutable y = a
            while y > b do
                bj_ii (float y) v
                saxpy 1.0 v (result.Ptr (offset + y * states)) v states
                for s = 0 to steps - 1 do // Integrate backwards over [y, y-1]
                    let t = float y - float s / float steps
                    // Hack: Fake limit from left
                    dV (if s = 0 then t - 1e-14 else t) v k1 //Tuple because of weird-ass parameter behaviour
                    sax h k1 k1 states
                    saxpy 0.5 k1 v tmp states
                    dV (t + h / 2.0) tmp k2
                    sax h k2 k2 states
                    saxpy 0.5 k2 v tmp states
                    dV (t + h / 2.0) tmp k3
                    sax h k3 k3 states
                    saxpy 1.0 k3 v tmp states
                    // Hack: Fake limit from right
                    dV (if s = steps - 1 then t + h + 1e-14 else t + h) tmp k4
                    sax h k4 k4 states
                    saxpy (1.0 / 6.0) k4 v tmp states
                    saxpy (1.0 / 3.0) k3 tmp tmp states
                    saxpy (1.0 / 3.0) k2 tmp tmp states
                    saxpy (1.0 / 6.0) k1 tmp v states
                arrayCopy v (result.Ptr (offset + y*states - states)) states
                y <- y - 1
        @> |> Compiler.DefineKernel 

    return Entry(fun program ->
        let worker = program.Worker
        let kernel = program.Apply kernel

        fun a b steps blockSize gridSize ->
            let n = (a + 1) * states * blockSize * gridSize
            use Va = worker.Malloc<float>(Array.zeroCreate states)
            use result = worker.Malloc(Array.zeroCreate n)
            
            let lp = LaunchParam (gridSize, blockSize)

            //worker.ProfilerStart()

            use start = worker.CreateEvent()
            use stop = worker.CreateEvent()
            worker.Synchronize()
            start.Record()
                      
            kernel.Launch lp a b steps Va.Ptr result.Ptr
            stop.Record()
            stop.Synchronize()

            //worker.ProfilerStop()
            let msec = Event.ElapsedMilliseconds(start, stop)
            let result = result.Gather()

            result, msec
        )
}