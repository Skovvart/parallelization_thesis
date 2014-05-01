module skovvart.actulus.cubase.TermInsurance

open skovvart.actulus.cubase.Common
open Alea.CUDA

type TermInsurance (planParams, paramCount) =
    inherit Plan<floatP>(50, 0, 1, planParams, paramCount) with
        [<ReflectedDefinition>]
        let b_0 t = conv 0
    
        [<ReflectedDefinition>]
        let mu_01 t age = GM t age
    
        [<ReflectedDefinition>]
        let bj_00 t = conv 0

        [<ReflectedDefinition>]
        let bj_01 t bdeath n = bdeath * indicator (t > conv 0) * indicator (t < n)

        override this.dV = <@ fun t (V:deviceptr<'T>) (planParams:deviceptr<'T>) (res:deviceptr<'T>) -> 
                let age = planParams.[0]
                let bdeath = planParams.[1]
                let n = planParams.[2]
                res.[0] <- (r t) * V.[0] - (b_0 t) - (mu_01 t age) * (conv 0 - V.[0] + (bj_01 t bdeath n))
            @>
        override this.bj_ii = <@ fun t (planParams:deviceptr<'T>) (res:deviceptr<'T>) -> 
                res.[0] <- bj_00 t
            @>

//let ti_bdeath = conv 1
//let ti_n = conv 35
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let ti_b_0 t = conv 0
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let ti_mu_01 t = GM t
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let ti_bj_00 t = conv 0
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let ti_bj_01 t = ti_bdeath * indicator (t > conv 0) * indicator (t < ti_n)
//
//let ti_dV = 
//    <@ fun (t:'T) (V:deviceptr<'T>) (res:deviceptr<'T>) -> 
//        res.[0] <- r t * V.[0] - ti_b_0 t - ti_mu_01 t * (conv 0 - V.[0] + ti_bj_01 t)
//    @>
//
//let ti_bj_ii =
//    <@ fun (t:'T) (res:deviceptr<'T>) -> 
//        res.[0] <- ti_bj_00 t
//    @>