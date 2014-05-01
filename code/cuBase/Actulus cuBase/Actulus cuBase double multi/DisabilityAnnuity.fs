module skovvart.actulus.cubase.DisabilityAnnuity

open skovvart.actulus.cubase.PlansCommon
open Alea.CUDA

let da_n = 35.0
let da_bdisabled = 1.0

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_b_0 t = 0.0

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_b_1 t = da_bdisabled * indicator (t > 0.0) * indicator (t < da_n)

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_GM01 t = 0.0006 + exp(ln10*(4.71609 - 10.0 + 0.06 * (age + t)))

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
let da_bj_00 t = 0.0

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_bj_01 t = 0.0

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_bj_02 t = 0.0

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_bj_11 t = 0.0

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let da_bj_12 t = 0.0


//Weird - only seems to accept up to two parameters, hence the Tuple
let da_dV = 
    <@ fun (t:float) (V:deviceptr<float>) (res:deviceptr<float>) -> 
        res.[0] <- r t * V.[0] - da_b_0 t - da_mu_01 t * (V.[1] - V.[0] + da_bj_01 t) - da_mu_02 t * (0.0 - V.[0] + da_bj_02 t)
        res.[1] <- r t * V.[1] - da_b_1 t - da_mu_12 t * (0.0 - V.[1] + da_bj_12 t)
    @>

let da_bj_ii =
    <@ fun (t:float) (res:deviceptr<float>) -> 
        res.[0] <- da_bj_00 t
        res.[1] <- da_bj_11 t
    @>