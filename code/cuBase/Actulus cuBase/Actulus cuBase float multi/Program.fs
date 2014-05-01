module skovvart.actulus.cubase.Main
open Alea.CUDA
open Alea.CUDA.Utilities
open skovvart.actulus.cubase.Common
open skovvart.actulus.cubase.PureEndowment
open skovvart.actulus.cubase.DeferredTemporaryLifeAnnuity
open skovvart.actulus.cubase.TemporaryLifeAnnuityPremium
open skovvart.actulus.cubase.TermInsurance
open skovvart.actulus.cubase.DisabilityAnnuity
open skovvart.actulus.cubase.DisabilityTermEnsurance
open skovvart.actulus.cubase.Device

let printResults states blocks threads (results:'T[]) = 
    let years = results.Length / states / blocks / threads
    let brange = blocks
    let trange = threads
    for b = 0 to brange - 1 do
        for t = 0 to trange - 1 do
            let offset = (t + b * threads) * states * years
            for i = 0 to years - 1 do
                printf "%3d:" i
                for s = 0 to states - 1 do
                    printf "  %20.16f" results.[offset + i*states + s]
                printfn ""
            printfn ""

let compile (p:Plan<_>) =
//    printfn "Compiling: %s" (p.GetType().Name)
    RK4_n p |> Compiler.load Worker.Default

//todo add variability
let genParams blocks threads planParams = 
    let amount = blocks*threads
    let variability p = 
        p |> Array.iteri (fun i -> fun t -> 
            let ageIndex = List.length planParams
            if (i % ageIndex) = 0 then p.[i] <- p.[i] - conv i/conv ageIndex)
        p
    let rec conc elems result times =
        match times with
        | 0 -> result
        | t -> conc elems (elems @ result) (times-1)
    conc planParams [] amount |> List.toArray |> variability

//untimed parameter gen
let default_params = [30; 1; 35]
let default_params_10 = default_params @ [10]
let d_default_params blocks threads = default_params |> List.map conv |> genParams blocks threads
let d_default_params_10 blocks threads = default_params_10 |> List.map conv |> genParams blocks threads

let PureEndowment blocks threads compile =
    let p = new PureEndowment(d_default_params blocks threads, default_params.Length)
    compile p, p.a, p.b, p.states

let DeferredTemporaryLifeAnnuity blocks threads compile =
    let p = new DeferredTemporaryLifeAnnuity(d_default_params_10 blocks threads, default_params_10.Length)
    compile p, p.a, p.b, p.states

let TemporaryLifeAnnuityPremium blocks threads compile =
    let p = new TemporaryLifeAnnuityPremium(d_default_params blocks threads, default_params.Length)
    compile p, p.a, p.b, p.states

let TermInsurance blocks threads compile =
    let p = new TermInsurance(d_default_params blocks threads, default_params.Length)
    compile p, p.a, p.b, p.states

let DisabilityAnnuity blocks threads compile =
    let p = new DisabilityAnnuity(d_default_params blocks threads, default_params.Length)
    compile p, p.a, p.b, p.states

let DisabilityTermInsurance blocks threads compile =
    let p = new DisabilityTermInsurance(d_default_params blocks threads, default_params.Length)
    compile p, p.a, p.b, p.states

let Plans blocks threads (compile:Plan<_>->'a) =
    [ PureEndowment blocks threads compile;
      DeferredTemporaryLifeAnnuity blocks threads compile;
      TemporaryLifeAnnuityPremium blocks threads compile;
      TermInsurance blocks threads compile;
      DisabilityAnnuity blocks threads compile;
      DisabilityTermInsurance blocks threads compile;
    ]

let run (program:Program<int -> int -> int -> int -> int -> 'T[] * float>) a b blocks threads =
    program.Run a b steps blocks threads

let unzip4 list =
    let l1 = list |> List.map (fun e -> 
        let (first, _, _, _) = e
        first)
    let l2 = list |> List.map (fun e -> 
        let (_, second, _, _) = e
        second)
    let l3 = list |> List.map (fun e -> 
        let (_, _, third, _ ) = e
        third)
    let l4 = list |> List.map (fun e -> 
        let (_, _, _, fourth) = e
        fourth)
    l1, l2, l3, l4

let runKernel plan blocks threads = 
    let kernel, a, b, states = plan
    let result, kernelTime = run kernel a b blocks threads
    kernelTime, result, states

let timeRunKernel plan blocks threads = 
    let execTime, (kernelTime, result, states) = time false "" (fun () -> runKernel plan blocks threads)
    execTime, kernelTime, result, states

let runAll plans print blocks threads =
    let execTimes, kernelTimes, results, states = List.map (fun plan -> timeRunKernel plan blocks threads) plans |> unzip4
    if print then
        List.iter2 (fun r -> fun s -> printResults s blocks threads r) results states
    kernelTimes |> List.sum, execTimes |> List.sum
    

//Force cuBase to do it's evaluation-prompt before timing begins
let untimed_cuBaseEvaluationPrompt = 
    let activePrompt = cuda { 
        let! kernel = <@ fun (empty:int) -> () @> |> Compiler.DefineKernel
        return Entry(fun program -> 
            let kernel = program.Apply kernel 
            ())
    }
    activePrompt |> Compiler.load Worker.Default |> ignore

let executeAverage blocks threads iterations print verbose = 
    time verbose "Program" (fun() -> 
        if verbose then
            printfn "Launching kernels with %d blocks and %d threads for a total of %d calculations, %d iterations" blocks threads (blocks*threads) iterations
        let compTime, plans = time verbose "Kernel compilation" (fun () -> Plans blocks threads compile)
        let mutable kernelTime = 0.0
        let mutable execTime = 0.0
        for i = 0 to iterations-1 do
            let kTime, eTime = runAll plans (print i) blocks threads
            kernelTime <- kernelTime + kTime
            execTime <- execTime + eTime
        let avgKernelTime = (kernelTime/float iterations)
        let avgExecTime = (execTime/float iterations)
        let memCpyTime = float execTime - kernelTime
        let avgMemCpyTime = (memCpyTime/float iterations)
        if verbose then
            printfn "Average sum of kernel-times: %.3fms" avgKernelTime
            printfn "Average execution-time: %.3fms" avgExecTime
            printfn "Average memory copying time: %.3f ms" avgMemCpyTime
        avgKernelTime
    )

let testPerformance blocks threads iterations print verbose = 
    let results = new System.Collections.Generic.Dictionary<int, System.Collections.Generic.Dictionary<int, float>>()
    for b in blocks do
        results.Add(b, new System.Collections.Generic.Dictionary<int, float>())
        for t in threads do
            try
                let _, avgKernelTime = executeAverage b t iterations print verbose
                results.[b].Add (t, avgKernelTime)
            with
            | :? System.Exception -> printfn "Out of memory (Blocks:%d, Threads:%d)\n" b t
    
    for b in results.Keys do
        for t in results.[b].Keys do
            let kernelRuntime = results.[b].[t]
            let calcsPerMs = float (b*t)/kernelRuntime
            printfn "(%d, %d): %f; %.3f" b t kernelRuntime calcsPerMs

let testResults blocks threads =
    let plans = Plans blocks threads compile
    try
        runAll plans true blocks threads |> ignore
    with
    | :? System.Exception -> printfn "Out of memory (Blocks:%d, Threads:%d)\n" blocks threads

[<EntryPoint>]
let main argv = 
    let print i = false//i=0
    let verbose = false
    let iterations = 100
    let blockSizes = [1; 14; 14*5; 14*10; 14*20; 14*25; 14*30]
    let threadSizes = [1; 8; 16; 32; 64; 128; 256; 512; 1024]

    if verbose then
        printDeviceInfo()

    time true "Program" (fun () -> testPerformance blockSizes threadSizes iterations print verbose) |> ignore
    //testResults 14 1
    0