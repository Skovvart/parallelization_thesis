module skovvart.actulus.cubase.TemporaryLifeAnnuityPremium

open skovvart.actulus.cubase.Common
open Alea.CUDA

type TemporaryLifeAnnuityPremium (planParams, paramCount) =
    inherit Plan<floatP>(50, 0, 1, planParams, paramCount) with
        [<ReflectedDefinition>]
        let b_0 t bpremium n = -bpremium * indicator (t >= conv 0) * indicator (t < n)
    
        [<ReflectedDefinition>]
        let mu_01 t age = GM t age
    
        [<ReflectedDefinition>]
        let bj_00 t = conv 0

        [<ReflectedDefinition>]
        let bj_01 t = conv 0

        override this.dV = <@ fun t (V:deviceptr<'T>) (planParams:deviceptr<'T>) (res:deviceptr<'T>) -> 
                let age = planParams.[0]
                let bpremium = planParams.[1]
                let n = planParams.[2]
                res.[0] <- (r t) * V.[0] - (b_0 t bpremium n) - (mu_01 t age) * (conv 0 - V.[0] + (bj_01 t))
            @>
        override this.bj_ii = <@ fun t (planParams:deviceptr<'T>) (res:deviceptr<'T>) -> 
                res.[0] <- bj_00 t
            @>

//let tlap_bpremium = conv 1
//let tlap_n = conv 35

//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let tlap_b_0 t = -tlap_bpremium * indicator (t >= conv 0) * indicator (t < tlap_n)
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let tlap_mu_01 t = GM t
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
//let tlap_bj_00 t = conv 0
//
//[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations (or maybe just in cuBase?) @>
//let tlap_bj_01 t = conv 0
//
//let tlap_dV = 
//    <@ fun (t:'T) (V:deviceptr<'T>) (res:deviceptr<'T>) -> 
//        res.[0] <- (r t) * V.[0] - (tlap_b_0 t) - (tlap_mu_01 t) * (conv 0 - V.[0] + (tlap_bj_01 t))
//    @>
//
//let tlap_bj_ii =
//    <@ fun (t:'T) (res:deviceptr<'T>) -> 
//        res.[0] <- tlap_bj_00 t
//    @>