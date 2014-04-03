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

//let printResults states results = 
//    Array.iteri (fun i -> fun e -> 
//        printfn "(%2d, %d): %20.16f" (i/states) (i%states) e) results
//    printfn ""

let stdFormatPrint states (results:float[]) = 
    let years = results.Length / states
    for i = 0 to years - 1 do
        printf "%3d:" i
        for s = 0 to states - 1 do
            printf "  %20.16f" results.[i*states + s]
        printfn ""
    printfn ""

let runTest dV bj_ii a b states =
    use program = RK4_n dV bj_ii states |> Compiler.load Worker.Default
    program.Run a b steps

let runPureEndowment =
    runTest pe_dV pe_bj_ii 40 0  1

let runDeferredTemporaryLifeAnnuity =
    runTest dtla_dV dtla_bj_ii 50 0 1

let runTemporaryLifeAnnuityPremium =
    runTest tlap_dV tlap_bj_ii 50 0 1

let runTermInsurance =
    runTest ti_dV ti_bj_ii 50 0 1

let runDisabilityAnnuity =
    runTest da_dV da_bj_ii 50 0 2

let runDisabilityTermInsurance =
    runTest dti_dV dti_bj_ii 50 0 2

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
    printfn "GPU: %s" Worker.Default.Device.Name
    let mutable totalTime = 0.0
    let iterations = 1
    for i = 0 to iterations-1 do
        let time = runAll true
        totalTime <- totalTime + time

    printfn "Tests finished in %fms" (totalTime/float iterations)
    0