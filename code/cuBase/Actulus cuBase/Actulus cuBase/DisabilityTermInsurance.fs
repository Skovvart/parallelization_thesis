module skovvart.actulus.cubase.DisabilityTermEnsurance

open skovvart.actulus.cubase.PlansCommon
open Alea.CUDA

let dti_n = 35.0
let dti_bdisabled = 1.0

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dti_b_0 t = 0.0

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dti_b_1 t = 0.0

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dti_GM01 t = 0.0006 + exp(ln10*(4.71609 - 10.0 + 0.06 * (age + t)))

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dti_GM02 t = GM t

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dti_GM12 t = GM t

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dti_mu_01 t = dti_GM01 t

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dti_mu_02 t = dti_GM02 t

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dti_mu_12 t = dti_GM12 t

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dti_bj_00 t = 0.0

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dti_bj_01 t = dti_bdisabled * indicator (t > 0.0) * indicator (t < dti_n)

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dti_bj_02 t = 0.0

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dti_bj_11 t = 0.0

[<ReflectedDefinition>] //Required to use methods (not constants) in <@ quotations @>
let dti_bj_12 t = 0.0


//Weird - only seems to accept up to two parameters, hence the Tuple
let dti_dV = 
    <@ fun (t:float) -> fun (V:deviceptr<float>, res:deviceptr<float>) -> 
        res.[0] <- r t * V.[0] - dti_b_0 t - dti_mu_01 t * (V.[1] - V.[0] + dti_bj_01 t) - dti_mu_02 t * (0.0 - V.[0] + dti_bj_02 t)
        res.[1] <- r t * V.[1] - dti_b_1 t - dti_mu_12 t * (0.0 - V.[1] + dti_bj_12 t)
    @>

let dti_bj_ii =
    <@ fun (t:float) -> fun (res:deviceptr<float>) -> 
        res.[0] <- dti_bj_00 t
        res.[1] <- dti_bj_11 t
    @>