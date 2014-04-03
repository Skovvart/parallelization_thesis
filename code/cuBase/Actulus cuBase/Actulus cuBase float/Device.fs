module skovvart.actulus.cubase.Device

open Microsoft.FSharp.Quotations
open Alea.CUDA
open Alea.CUDA.Utilities

[<ReflectedDefinition>]
let arrayCopy (src:deviceptr<float32>) (destination:deviceptr<float32>) length =
    for i = 0 to length - 1 do
        destination.[i] <- src.[i]

[<ReflectedDefinition>]
let sax a (x:deviceptr<float32>) (res:deviceptr<float32>) length =
    for i = 0 to length - 1 do
        res.[i] <- a * x.[i]

[<ReflectedDefinition>]
let saxpy a (x:deviceptr<float32>) (y:deviceptr<float32>) (res:deviceptr<float32>) length =
    for i = 0 to length - 1 do
        res.[i] <- a * x.[i] + y.[i]

let RK4_n dV bj_ii states = cuda {
    let! kernel =
        <@ fun (a:int) b steps (Va:deviceptr<float32>) (result:deviceptr<float32>) ->
            let offset = (blockIdx.x * blockDim.x + threadIdx.x) * states * (a + 1)
            let dV = %dV
            let bj_ii = %bj_ii

            let h   = -1.0f / float32 steps
            let k1  = __local__.Array(states) |> __array_to_ptr
            let k2  = __local__.Array(states) |> __array_to_ptr
            let k3  = __local__.Array(states) |> __array_to_ptr
            let k4  = __local__.Array(states) |> __array_to_ptr
            let tmp = __local__.Array(states) |> __array_to_ptr
            let v   = __local__.Array(states) |> __array_to_ptr
            
            arrayCopy Va (result.Ptr (offset + a * states)) states
            let mutable y = a
            while y > b do
                bj_ii (float32 y) v
                saxpy 1.0f v (result.Ptr (offset + y * states)) v states
                for s = 0 to steps - 1 do // Integrate backwards over [y, y-1]
                    let t = float32 y - float32 s / float32 steps
                    // Hack: Fake limit from left
                    dV (match s with | 0 -> t - 1e-5f | _ -> t) (v, k1) //Tuple because of weird-ass parameter behaviour
                    sax h k1 k1 states
                    saxpy 0.5f k1 v tmp states
                    dV (t + h / 2.0f) (tmp, k2)
                    sax h k2 k2 states
                    saxpy 0.5f k2 v tmp states
                    dV (t + h / 2.0f) (tmp, k3)
                    sax h k3 k3 states
                    saxpy 1.0f k3 v tmp states
                    // Hack: Fake limit from right
                    dV (match s with | _ when s = steps - 1 -> t + h + 1e-5f | _ -> t + h) (tmp, k4)
                    sax h k4 k4 states
                    saxpy (1.0f / 6.0f) k4 v tmp states
                    saxpy (2.0f / 6.0f) k3 tmp tmp states
                    saxpy (2.0f / 6.0f) k2 tmp tmp states
                    saxpy (1.0f / 6.0f) k1 tmp v states
                arrayCopy v (result.Ptr (offset + y*states - states)) states
                y <- y - 1
        @> |> Compiler.DefineKernel 

    return Entry(fun program ->
        let worker = program.Worker
        let kernel = program.Apply kernel

        fun a b steps ->
            let n = (a + 1) * states
            use Va = worker.Malloc(Array.zeroCreate states)
            use result = worker.Malloc(Array.zeroCreate n)

            let blockSize = 128
            let numSm = worker.Device.Attributes.MULTIPROCESSOR_COUNT
            // We tend to partition data so that each SM could handle
            // 16 blocks to hide the memory latency.
            // For more detail, please reference "SM Occupancy".
            //divup = divide and round up / ceiling
            let gridSize = min (numSm * 16) (divup n blockSize)
            // Now we know the launch shape, could create launching parameter.
            //printfn "Launching kernel with blocksize %d and gridsize %d" blockSize gridSize
            let lp = LaunchParam (gridSize, blockSize)

            use start = worker.CreateEvent()
            use stop = worker.CreateEvent()
            worker.Synchronize()
            start.Record()

            kernel.Launch lp a b steps Va.Ptr result.Ptr

            stop.Record()
            stop.Synchronize()
            let msec = Event.ElapsedMilliseconds(start, stop)

            let result = result.Gather()

            result, msec
        )
}