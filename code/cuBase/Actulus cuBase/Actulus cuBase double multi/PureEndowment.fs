module skovvart.actulus.cubase.PureEndowment

open skovvart.actulus.cubase.PlansCommon
open Alea.CUDA

[<ReflectedDefinition>] //Required to use in <@ quotations @>
let pe_b_0 t = 0.0

[<ReflectedDefinition>]
let pe_mu_01 t = GM t

[<ReflectedDefinition>]
let pe_bj_00 t = if t = pensiontime then bpension else 0.0

[<ReflectedDefinition>]
let pe_bj_01 t = 0.0

[<ReflectedDefinition>] //Required to use in <@ quotations @>
let pe_dV = 
    <@ fun (t:float) (V:deviceptr<float>) (res:deviceptr<float>) -> 
        res.[0] <- (r t) * V.[0] - (pe_b_0 t) - (pe_mu_01 t) * (0.0 - V.[0] + (pe_bj_01 t))
    @>

let pe_bj_ii =
    <@ fun t (res:deviceptr<float>) -> 
        res.[0] <- if t = pensiontime then bpension else pe_bj_00 t
    @>