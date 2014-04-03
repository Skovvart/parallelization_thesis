module skovvart.actulus.cubase.PlansCommon

let age = 30.0
let interestrate = 0.05
let bpension = 1.0
let pensiontime = 35.0

[<ReflectedDefinition>] //Required to use in <@ quotations @>
let indicator b = if b then 1.0 else 0.0

[<ReflectedDefinition>]
let r t = 0.05

let ln10 = 2.3025850929940459

[<ReflectedDefinition>]
let GM t = 0.0005 + exp (ln10 * (5.728 - 10.0 + 0.038*(age + t)))