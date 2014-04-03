﻿module skovvart.actulus.cubase.PureEndowment

open skovvart.actulus.cubase.PlansCommon
open Alea.CUDA

[<ReflectedDefinition>] //Required to use in <@ quotations @>
let pe_b_0 t = 0.0f

[<ReflectedDefinition>]
let pe_mu_01 t = GM t

[<ReflectedDefinition>]
let pe_bj_00 t = if t = pensiontime then bpension else 0.0f

[<ReflectedDefinition>]
let pe_bj_01 t = 0.0f

//Weird - only seems to accept up to two parameters, hence the Tuple
[<ReflectedDefinition>] //Required to use in <@ quotations @>
let pe_dV = 
    <@ fun (t:float32) -> fun (V:deviceptr<float32>, res:deviceptr<float32>) -> 
        res.[0] <- (r t) * V.[0] - (pe_b_0 t) - (pe_mu_01 t) * (0.0f - V.[0] + (pe_bj_01 t))
    @>

let pe_bj_ii =
    <@ fun t -> fun (res:deviceptr<float32>) -> 
        res.[0] <- if t = pensiontime then bpension else pe_bj_00 t
    @>