module skovvart.actulus.cubase.TemporaryLifeAnnuityPremium

open skovvart.actulus.cubase.PlansCommon
open Alea.CUDA

let tlap_n = 35.0f
let tlap_bpremium = 1.0f

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let tlap_b_0 t = -tlap_bpremium * indicator (t >= 0.0f) * indicator (t < tlap_n)

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let tlap_mu_01 t = GM t

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let tlap_bj_00 t = 0.0f

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations (or maybe just in cuBase?) @>
let tlap_bj_01 t = 0.0f


//Weird - only seems to accept up to two parameters, hence the Tuple
let tlap_dV = 
    <@ fun (t:float32) -> fun (V:deviceptr<float32>, res:deviceptr<float32>) -> 
        res.[0] <- (r t) * V.[0] - (tlap_b_0 t) - (tlap_mu_01 t) * (0.0f - V.[0] + (tlap_bj_01 t))
    @>

let tlap_bj_ii =
    <@ fun (t:float32) -> fun (res:deviceptr<float32>) -> 
        res.[0] <- tlap_bj_00 t
    @>