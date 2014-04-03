module skovvart.actulus.cubase.TermInsurance

open skovvart.actulus.cubase.PlansCommon
open Alea.CUDA

let ti_n = 35.0f
let ti_bdeath = 1.0f

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let ti_b_0 t = 0.0f

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let ti_mu_01 t = GM t

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let ti_bj_00 t = 0.0f

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let ti_bj_01 t = ti_bdeath * indicator (t > 0.0f) * indicator (t < ti_n)


//Weird - only seems to accept up to two parameters, hence the Tuple
let ti_dV = 
    <@ fun (t:float32) -> fun (V:deviceptr<float32>, res:deviceptr<float32>) -> 
        res.[0] <- r t * V.[0] - ti_b_0 t - ti_mu_01 t * (0.0f - V.[0] + ti_bj_01 t)
    @>

let ti_bj_ii =
    <@ fun (t:float32) -> fun (res:deviceptr<float32>) -> 
        res.[0] <- ti_bj_00 t
    @>