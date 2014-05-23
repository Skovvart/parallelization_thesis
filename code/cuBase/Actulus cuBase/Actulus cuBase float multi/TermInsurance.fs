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