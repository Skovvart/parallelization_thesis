module skovvart.actulus.cubase.DeferredTemporaryLifeAnnuity

open skovvart.actulus.cubase.PlansCommon
open Alea.CUDA

let dtla_m = 35.0f
let dtla_n = 10.0f

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dtla_b_0 t = bpension * indicator (t > dtla_m) * indicator (t < dtla_m + dtla_n)

[<ReflectedDefinition>] //Required to use in <@ quotations @>
let dtla_mu_01 t = GM t

[<ReflectedDefinition>] //Required to use in <@ quotations @>
let dtla_bj_00 t = 0.0f

[<ReflectedDefinition>] //Required to use in <@ quotations @>
let dtla_bj_01 t = 0.0f

//Weird - only seems to accept up to two parameters, hence the Tuple
let dtla_dV = 
    <@ fun (t:float32) -> fun (V:deviceptr<float32>, res:deviceptr<float32>) -> 
        res.[0] <- (r t) * V.[0] - (dtla_b_0 t) - (dtla_mu_01 t) * (0.0f - V.[0] + (dtla_bj_01 t))
    @>

let dtla_bj_ii =
    <@ fun (t:float32) -> fun (res:deviceptr<float32>) -> 
        res.[0] <- dtla_bj_00 t
    @>