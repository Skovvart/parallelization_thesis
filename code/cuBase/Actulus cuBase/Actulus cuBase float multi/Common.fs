module skovvart.actulus.cubase.Common

open Alea.CUDA

let steps = 100

//Plan "interface"
[<AbstractClass>]
type Plan<'T>(a:int, b:int, states:int, planParams:'T[], paramCount:int) =
    member this.a = a
    member this.b = b
    member this.states = states
    member this.planParams = planParams
    member this.paramCount = paramCount
    abstract member dV : Quotations.Expr<('T -> deviceptr<'T> -> deviceptr<'T> -> deviceptr<'T> -> unit)>
    abstract member bj_ii : Quotations.Expr<('T -> deviceptr<'T> -> deviceptr<'T> -> unit)>

//type converters - update both floatP and conv
type floatP = float
[<ReflectedDefinition>]
let inline conv v  = 
    float v

let intConv = 
    <@ conv @>

let floatConv = //gimmicky way to avoid having to manually change this when changing between doubles/floats
    let convResType = (conv 0).GetType()
    if convResType = typeof<float> then
        <@ (fun f -> if f = 1e-5 then conv 1e-14 else conv f) @>
    else
        <@ conv @>

//common kernel methods and constants
[<ReflectedDefinition>]
let indicator b = if b then conv 1 else conv 0

[<ReflectedDefinition>]
let r t = conv 0.05

let ln10 = conv 2.3025850929940459

[<ReflectedDefinition>]
let GM t age = conv 0.0005 + exp (ln10 * (conv 5.728 - (conv 10) + conv 0.038*(age + t)))

//Util
let time print desc prog = 
    if print then
        printfn "%s started" desc
    let worker = Worker.Default
    use start = worker.CreateEvent()
    use stop = worker.CreateEvent()
    worker.Synchronize()
    start.Record()
    let res = prog()
    stop.Record()
    stop.Synchronize()
    let msec = Event.ElapsedMilliseconds(start, stop)
    if print then
        printfn "%s finished in %.3f ms" desc msec
        printfn ""
    msec, res

let printDeviceInfo () = 
    let d = Worker.Default.Device
    let a = d.Attributes
    
    printfn "GPU: %s - Compute version %d.%d (mode: %d), %.3fGHz (memory: %.3fGHz)" d.Name a.COMPUTE_CAPABILITY_MAJOR a.COMPUTE_CAPABILITY_MINOR a.COMPUTE_MODE (float a.CLOCK_RATE/1e6) (float a.MEMORY_CLOCK_RATE/1e6)
    printfn "SMs: %d, Threads per Block/SM: %d/%d, Concurrent kernels: %d" a.MULTIPROCESSOR_COUNT a.MAX_THREADS_PER_BLOCK a.MAX_THREADS_PER_MULTIPROCESSOR a.CONCURRENT_KERNELS
    printfn "Block dims (%d,%d,%d)" a.MAX_BLOCK_DIM_X a.MAX_BLOCK_DIM_Y a.MAX_BLOCK_DIM_Z
    printfn "Thread dims (%d,%d,%d)" a.MAX_GRID_DIM_X a.MAX_GRID_DIM_Y a.MAX_GRID_DIM_Z
    printfn "Warp size: %d" a.WARP_SIZE
    printfn "Memory - Constant: %d, shared per block: %d, registers per block: %d, L2 cache size: %d" a.TOTAL_CONSTANT_MEMORY a.SHARED_MEMORY_PER_BLOCK a.REGISTERS_PER_BLOCK a.L2_CACHE_SIZE
    printfn ""