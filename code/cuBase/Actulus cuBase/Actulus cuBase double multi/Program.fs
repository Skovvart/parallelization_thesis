module skovvart.actulus.cubase.Main

open Alea.CUDA
open Alea.CUDA.Utilities
open skovvart.actulus.cubase.PureEndowment
open skovvart.actulus.cubase.DeferredTemporaryLifeAnnuity
open skovvart.actulus.cubase.TemporaryLifeAnnuityPremium
open skovvart.actulus.cubase.TermInsurance
open skovvart.actulus.cubase.DisabilityAnnuity
open skovvart.actulus.cubase.DisabilityTermEnsurance
open skovvart.actulus.cubase.Device

let steps = 100
let blocks = 64
let threads = 1536

let stdFormatPrint states (results:float[]) = 
    let years = results.Length / states / blocks / threads
    //blocks
    for b = 0 to 1 - 1 do
        //threads
        for t = 0 to 1 - 1 do
            let offset = (t + b * threads) * states * years
            for i = 0 to years - 1 do
                printf "%3d:" i
                for s = 0 to states - 1 do
                    printf "  %20.16f" results.[offset + i*states + s]
                printfn ""
            printfn ""

let run dV bj_ii a b states =
    use program = RK4_n dV bj_ii states |> Compiler.load Worker.Default
    program.Run a b steps blocks threads

let runPureEndowment =
    run pe_dV pe_bj_ii 40 0 1

let runDeferredTemporaryLifeAnnuity =
    run dtla_dV dtla_bj_ii 50 0 1

let runTemporaryLifeAnnuityPremium =
    run tlap_dV tlap_bj_ii 50 0 1

let runTermInsurance =
    run ti_dV ti_bj_ii 50 0 1

let runDisabilityAnnuity =
    run da_dV da_bj_ii 50 0 2

let runDisabilityTermInsurance =
    run dti_dV dti_bj_ii 50 0 2

let runAll print =
    let totalTime = 0.0
    let results, time = runPureEndowment
    if print then
        results |> stdFormatPrint 1
    let totalTime = totalTime + time

    let results, time = runDeferredTemporaryLifeAnnuity
    if print then
        results |> stdFormatPrint 1
    let totalTime = totalTime + time

    let results, time = runTemporaryLifeAnnuityPremium
    if print then
        results |> stdFormatPrint 1
    let totalTime = totalTime + time

    let results, time = runTermInsurance
    if print then
        results |> stdFormatPrint 1
    let totalTime = totalTime + time

    let results, time = runDisabilityAnnuity
    if print then
        results |> stdFormatPrint 2
    let totalTime = totalTime + time

    let results, time = runDisabilityTermInsurance
    if print then
        results |> stdFormatPrint 2
    totalTime + time

[<EntryPoint>]
let main argv = 
    printfn "GPU: %s, SMs: %d, Threads per SM: %d" Worker.Default.Device.Name Worker.Default.Device.Attributes.MULTIPROCESSOR_COUNT Worker.Default.Device.Attributes.MAX_THREADS_PER_MULTIPROCESSOR
    let mutable totalTime = 0.0
    let iterations = 1000
    for i = 0 to iterations-1 do
        let time = runAll (i=0)
        totalTime <- totalTime + time

    printfn "Tests finished in %fms" (totalTime/float iterations)
    0