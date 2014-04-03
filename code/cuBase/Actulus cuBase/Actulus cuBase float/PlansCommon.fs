module skovvart.actulus.cubase.PlansCommon

let age = 30.0f
let interestrate = 0.05f
let bpension = 1.0f
let pensiontime = 35.0f

[<ReflectedDefinition>] //Required to use in <@ quotations @>
let indicator b = if b then 1.0f else 0.0f

[<ReflectedDefinition>]
let r t = 0.05f

let ln10 = 2.3025850929940459f

[<ReflectedDefinition>]
let GM t = 0.0005f + exp (ln10 * (5.728f - 10.0f + 0.038f*(age + t)))