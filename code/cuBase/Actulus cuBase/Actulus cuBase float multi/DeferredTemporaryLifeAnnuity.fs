module skovvart.actulus.cubase.DeferredTemporaryLifeAnnuity

open skovvart.actulus.cubase.Common
open Alea.CUDA

type DeferredTemporaryLifeAnnuity (planParams, paramCount) =
    inherit Plan<floatP>(50, 0, 1, planParams, paramCount) with
        [<ReflectedDefinition>]
        let b_0 t pension m n = pension * indicator (t > m) * indicator (t < m + n)
        ////let dtla_b_0 t = dtla_bpension * indicator (t > dtla_m) * indicator (t < dtla_m + dtla_n)

        [<ReflectedDefinition>]
        let mu_01 t age = GM t age

        [<ReflectedDefinition>]
        let bj_00 t = conv 0

        [<ReflectedDefinition>]
        let bj_01 t = conv 0

        override this.dV = <@ fun t (V:deviceptr<'T>) (planParams:deviceptr<'T>) (res:deviceptr<'T>) -> 
                let age = planParams.[0]
                let pension = planParams.[1]
                let m = planParams.[2]
                let n = planParams.[3]
                res.[0] <- (r t) * V.[0] - (b_0 t pension m n) - (mu_01 t age) * (conv 0 - V.[0] + (bj_01 t))
            @>
        override this.bj_ii = <@ fun t (planParams:deviceptr<'T>) (res:deviceptr<'T>) -> 
                res.[0] <- bj_00 t
            @>


//let dtla_m = conv 35
//let dtla_n = conv 10
//let dtla_bpension = conv 1
//
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let dtla_b_0 t = dtla_bpension * indicator (t > dtla_m) * indicator (t < dtla_m + dtla_n)
//
//[<ReflectedDefinition>] //Required to use in <@ quotations @>
//let dtla_mu_01 t = GM t
//
//[<ReflectedDefinition>] //Required to use in <@ quotations @>
//let dtla_bj_00 t = conv 0
//
//[<ReflectedDefinition>] //Required to use in <@ quotations @>
//let dtla_bj_01 t = conv 0
//
//let dtla_dV = 
//    <@ fun (t:'T) (V:deviceptr<'T>) (res:deviceptr<'T>) -> 
//        res.[0] <- (r t) * V.[0] - (dtla_b_0 t) - (dtla_mu_01 t) * (conv 0 - V.[0] + (dtla_bj_01 t))
//    @>
//
//let dtla_bj_ii =
//    <@ fun (t:'T) (res:deviceptr<'T>) -> 
//        res.[0] <- dtla_bj_00 t
//    @>