module skovvart.actulus.cubase.DisabilityAnnuity

open skovvart.actulus.cubase.PlansCommon
open Alea.CUDA

let da_n = 35.0f
let da_bdisabled = 1.0f

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_b_0 t = 0.0f

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_b_1 t = da_bdisabled * indicator (t > 0.0f) * indicator (t < da_n)

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_GM01 t = 0.0006f + exp(ln10*(4.71609f - 10.0f + 0.06f * (age + t)))

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_GM02 t = GM t

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_GM12 t = GM t

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_mu_01 t = da_GM01 t

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_mu_02 t = da_GM02 t

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_mu_12 t = da_GM12 t

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_bj_00 t = 0.0f

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_bj_01 t = 0.0f

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_bj_02 t = 0.0f

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_bj_11 t = 0.0f

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_bj_12 t = 0.0f


//Weird - only seems to accept up to two parameters, hence the Tuple
let da_dV = 
    <@ fun (t:float32) -> fun (V:deviceptr<float32>, res:deviceptr<float32>) -> 
        res.[0] <- r t * V.[0] - da_b_0 t - da_mu_01 t * (V.[1] - V.[0] + da_bj_01 t) - da_mu_02 t * (0.0f - V.[0] + da_bj_02 t)
        res.[1] <- r t * V.[1] - da_b_1 t - da_mu_12 t * (0.0f - V.[1] + da_bj_12 t)
    @>

let da_bj_ii =
    <@ fun (t:float32) -> fun (res:deviceptr<float32>) -> 
        res.[0] <- da_bj_00 t
        res.[1] <- da_bj_11 t
    @>