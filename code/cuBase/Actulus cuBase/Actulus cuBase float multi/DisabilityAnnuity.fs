module skovvart.actulus.cubase.DisabilityAnnuity

open skovvart.actulus.cubase.Common
open Alea.CUDA

type DisabilityAnnuity (planParams, paramCount) =
    inherit Plan<floatP>(50, 0, 2, planParams, paramCount) with
        [<ReflectedDefinition>]
        let b_0 t = conv 0

        [<ReflectedDefinition>]
        let b_1 t bdisabled n = bdisabled * indicator (t > conv 0) * indicator (t < n)
        //planParams.[1] * indicator (t > conv 0) * indicator (t < planParams.[2])

        [<ReflectedDefinition>]
        let GM01 t age = conv 0.0006 + exp(ln10*(conv 4.71609 - conv 10 + conv 0.06 * (age + t)))

        [<ReflectedDefinition>]
        let GM02 t age = GM t age

        [<ReflectedDefinition>]
        let GM12 t age = GM t age

        [<ReflectedDefinition>]
        let mu_01 t age = GM01 t age

        [<ReflectedDefinition>]
        let mu_02 t age = GM02 t age

        [<ReflectedDefinition>]
        let mu_12 t age = GM12 t age

        [<ReflectedDefinition>]
        let bj_00 t = conv 0

        [<ReflectedDefinition>]
        let bj_01 t = conv 0

        [<ReflectedDefinition>]
        let bj_02 t = conv 0

        [<ReflectedDefinition>]
        let bj_11 t = conv 0

        [<ReflectedDefinition>]
        let bj_12 t = conv 0

        override this.dV = <@ fun (t:'T) (V:deviceptr<'T>) (planParams:deviceptr<'T>) (res:deviceptr<'T>) -> 
                let age = planParams.[0]
                let bdisabled = planParams.[1]
                let n = planParams.[2]
                res.[0] <- r t * V.[0] - b_0 t - mu_01 t age * (V.[1] - V.[0] + bj_01 t) - mu_02 t age * (conv 0 - V.[0] + bj_02 t)
                res.[1] <- r t * V.[1] - b_1 t bdisabled n - mu_12 t age * (conv 0 - V.[1] + bj_12 t)
            @>
        override this.bj_ii = <@ fun t (planParams:deviceptr<'T>) (res:deviceptr<'T>) -> 
                res.[0] <- bj_00 t
                res.[1] <- bj_11 t
            @>

//let da_bdisabled = conv 1
//let da_n = conv 35
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let da_b_0 t = conv 0
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let da_b_1 t = da_bdisabled * indicator (t > conv 0) * indicator (t < da_n)
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let da_GM01 t = conv 0.0006 + exp(ln10*(conv 4.71609 - conv 10 + conv 0.06 * (age + t)))
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let da_GM02 t = GM t
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let da_GM12 t = GM t
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let da_mu_01 t = da_GM01 t
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let da_mu_02 t = da_GM02 t
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let da_mu_12 t = da_GM12 t
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let da_bj_00 t = conv 0
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let da_bj_01 t = conv 0
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let da_bj_02 t = conv 0
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let da_bj_11 t = conv 0
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let da_bj_12 t = conv 0
//
//let da_dV = 
//    <@ fun (t:'T) (V:deviceptr<'T>) (res:deviceptr<'T>) -> 
//        res.[0] <- r t * V.[0] - da_b_0 t - da_mu_01 t * (V.[1] - V.[0] + da_bj_01 t) - da_mu_02 t * (conv 0 - V.[0] + da_bj_02 t)
//        res.[1] <- r t * V.[1] - da_b_1 t - da_mu_12 t * (conv 0 - V.[1] + da_bj_12 t)
//    @>
//
//let da_bj_ii =
//    <@ fun (t:'T) (res:deviceptr<'T>) -> 
//        res.[0] <- da_bj_00 t
//        res.[1] <- da_bj_11 t
//    @>