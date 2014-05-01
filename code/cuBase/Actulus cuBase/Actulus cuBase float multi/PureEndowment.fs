module skovvart.actulus.cubase.PureEndowment

open skovvart.actulus.cubase.Common
open Alea.CUDA

type PureEndowment (planParams, paramCount) =
    inherit Plan<floatP>(40, 0, 1, planParams, paramCount) with
        [<ReflectedDefinition>]
        let b_0 t = conv 0
    
        [<ReflectedDefinition>]
        let mu_01 t age = GM t age
    
        [<ReflectedDefinition>]
        let bj_00 t pensiontime bpension = if t = pensiontime then bpension else conv 0
        //let bj_00 t = if t = pensiontime then bpension else conv 0

        [<ReflectedDefinition>]
        let bj_01 t = conv 0

        override this.dV = <@ fun t (V:deviceptr<'T>) (planParams:deviceptr<'T>) (res:deviceptr<'T>) -> 
                let age = planParams.[0]
                res.[0] <- (r t) * V.[0] - (b_0 t) - (mu_01 t age) * (conv 0 - V.[0] + (bj_01 t))
            @>
        override this.bj_ii = <@ fun t (planParams:deviceptr<'T>) (res:deviceptr<'T>) -> 
                let pensiontime = planParams.[2]
                let bpension = planParams.[1]
                res.[0] <- bj_00 t pensiontime bpension
            @>
//
//[<ReflectedDefinition>] //Required to use in <@ quotations @>
//let pe_b_0 t = conv 0
//
//[<ReflectedDefinition>]
//let pe_mu_01 t = GM t
//
//[<ReflectedDefinition>]
//let pe_bj_00 t = if t = pensiontime then bpension else conv 0
//
//[<ReflectedDefinition>]
//let pe_bj_01 t = conv 0
//
//let pe_dV = 
//    <@ fun (t:'T) (V:deviceptr<'T>) (res:deviceptr<'T>) -> 
//        res.[0] <- (r t) * V.[0] - (pe_b_0 t) - (pe_mu_01 t) * (conv 0 - V.[0] + (pe_bj_01 t))
//    @>
//
//let pe_bj_ii =
//    <@ fun t (res:deviceptr<'T>) -> 
//        res.[0] <- pe_bj_00 t
//    @>