module skovvart.actulus.cubase.Host

open skovvart.actulus.cubase.PlansCommon
open skovvart.actulus.cubase.PureEndowment
open skovvart.actulus.cubase.DeferredTemporaryLifeAnnuity
open skovvart.actulus.cubase.TemporaryLifeAnnuityPremium
open skovvart.actulus.cubase.TermInsurance
open skovvart.actulus.cubase.DisabilityAnnuity
open skovvart.actulus.cubase.DisabilityTermEnsurance

let flatten (A:'a[,]) = A |> Seq.cast<'a>

let getRow r (A:_[,]) = flatten A.[r..r,*] |> Seq.toArray

let sax a x = Array.map (fun e -> a * e) x

let saxpy a x y = Array.map2 (fun xe -> fun ye -> a * xe + ye) x y

let RK4_n dV bj_ii a b steps (Va:float[]) =
    let n = Va.Length
    let mutable result = Array2D.zeroCreate (a - b + 1) n
    let h = -1.0 / (float steps)
    for i = 0 to n - 1 do
        result.[a - b, i] <- Va.[i]
    let mutable k1 = Array.zeroCreate n
    let mutable k2 = Array.zeroCreate n
    let mutable k3 = Array.zeroCreate n
    let mutable k4 = Array.zeroCreate n
    let mutable tmp = Array.zeroCreate n
    let mutable v = Array.zeroCreate n
    let mutable y = a
    while y > b do
        bj_ii (float y) v
        let saxpyres = saxpy 1.0 v (getRow (y - b) result)
        for i = 0 to n - 1 do
            v.[i] <- saxpyres.[i]
        for s = 0 to steps - 1 do // Integrate backwards over [y, y-1]
            let t = (float y) - (float s) / (float steps)
            // Hack: Fake limit from left
            dV (match s with | 0 -> t - 1e-14 | _ -> t) v k1
            k1 <- sax h k1
            tmp <- saxpy 0.5 k1 v
            dV (t + h / 2.0) tmp k2
            k2 <- sax h k2
            tmp <- saxpy 0.5 k2 v
            dV (t + h / 2.0) tmp k3
            k3 <- sax h k3
            tmp <- saxpy 1.0 k3 v
            // Hack: Fake limit from right
            dV (match s with | _ when s = steps - 1 -> t + h + 1e-14 | _ -> t + h) tmp k4
            v <- saxpy (1.0 / 6.0) k1 (saxpy (2.0 / 6.0) k2 (saxpy (2.0 / 6.0) k3 (saxpy (1.0 / 6.0) (sax h k4) v)))
        for i = 0 to n - 1 do
            result.[y - 1 - b, i] <- v.[i]
        y <- y - 1
    result

let indicator bool = if bool then 1.0 else 0.0
let ln10 = 2.3025850929940459
let GM t = 0.0005 + exp(ln10*(5.728 - 10.0 + 0.038 * (age + t)))
let r t = interestrate

let rFsa t = 
    let ts = [| 0.25; 0.5; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; 9.0; 10.0; 11.0; 12.0; 13.0; 14.0; 15.0; 16.0; 17.0; 18.0; 19.0; 20.0; 21.0; 22.0; 23.0; 24.0; 25.0; 26.0; 27.0; 28.0; 29.0; 30.0 |]
    let rs = [| 1.146677033; 1.146677033; 1.146677033; 1.340669678; 1.571952911; 1.803236144; 2.034519377; 2.265802610; 2.497085843; 2.584085843; 2.710085843; 2.805085843; 2.871485843; 2.937885843; 3.004285843; 3.070685843; 3.137085843; 3.136485843; 3.135885843; 3.135285843; 3.134685843; 3.134085843; 3.113185843; 3.092285843; 3.071385843; 3.050485843; 3.029585843; 3.008685843; 2.987785843; 2.966885843; 2.945985843; 2.925085843 |]
    // Requires ts non-empty and elements strictly increasing.
    let last = ts.Length - 1
    if t <= ts.[0] then log (1.0 + rs.[0] / 100.0) else
        if t >= ts.[last] then log (1.0 + rs.[last] / 100.0) else
            let mutable a = 0
            let mutable b = last
            // Now a < b (bcs. ts must have more than 1 element) and ts[a] < t < ts[b]
            while a + 1 < b do
                // Now a < b and ts[a] <= t < ts[b]
                let i = (a + b) / 2
                if ts.[i] <= t then
                    a <- i
                else // t < ts[i]
                    b <- i
            // Now a+1>=b and ts[a] <= t < ts[b]; so a!=b and hence a+1 == b <= last
            let m = a
            let tm = ts.[m]
            let tm1 = ts.[m + 1]
            let rm = rs.[m] / 100.0
            let rm1 = rs.[m + 1] / 100.0
            let Rt = (rm * (tm1 - t) + rm1 * (t - tm)) / (tm1 - tm)
            log(1.0 + Rt) + t / (tm1 - tm) * (rm1 - rm) / (1.0 + Rt)

let PureEndowment = 
    RK4_n 
        (fun t -> fun V -> fun res -> res.[0] <- (r t) * V.[0] - (pe_b_0 t) - (pe_mu_01 t) * (0.0 - V.[0] + (pe_bj_01 t))) 
        (fun t -> fun res -> res.[0] <- pe_bj_00 t) 
        40 0 100 (Array.zeroCreate 1)

let DeferredTemporaryLifeAnnuity = 
    RK4_n 
        (fun t -> fun V -> fun res -> res.[0] <- (r t) * V.[0] - (dtla_b_0 t) - (dtla_mu_01 t) * (0.0 - V.[0] + (dtla_bj_01 t))) 
        (fun t -> fun res -> res.[0] <- dtla_bj_00 t) 
        50 0 100 (Array.zeroCreate 1)

let TemporaryLifeAnnuityPremium = 
    RK4_n 
        (fun t -> fun V -> fun res -> res.[0] <- (r t) * V.[0] - (tlap_b_0 t) - (tlap_mu_01 t) * (0.0 - V.[0] + (tlap_bj_01 t))) 
        (fun t -> fun res -> res.[0] <- tlap_bj_00 t) 
        50 0 100 (Array.zeroCreate 1)

let TermInsurance = 
    RK4_n 
        (fun t -> fun V -> fun res -> res.[0] <- r t * V.[0] - ti_b_0 t - ti_mu_01 t * (0.0 - V.[0] + ti_bj_01 t)) 
        (fun t -> fun res -> res.[0] <- ti_bj_00 t) 
        50 0 100 (Array.zeroCreate 1)

let DisabilityAnnuity = 
    RK4_n 
        (fun t -> fun V -> fun res -> 
                res.[0] <- r t * V.[0] - da_b_0 t - da_mu_01 t * (V.[1] - V.[0] + da_bj_01 t) - da_mu_02 t * (0.0 - V.[0] + da_bj_02 t)
                res.[1] <- r t * V.[1] - da_b_1 t - da_mu_12 t * (0.0 - V.[1] + da_bj_12 t))
        (fun t -> fun res -> 
            res.[0] <- da_bj_00 t
            res.[1] <- da_bj_11 t)
        50 0 100 (Array.zeroCreate 2)

let DisabilityTermInsurance = 
    RK4_n 
        (fun t -> fun V -> fun res -> 
            res.[0] <- r t * V.[0] - dti_b_0 t - dti_mu_01 t * (V.[1] - V.[0] + dti_bj_01 t) - dti_mu_02 t * (0.0 - V.[0] + dti_bj_02 t)
            res.[1] <- r t * V.[1] - dti_b_1 t - dti_mu_12 t * (0.0 - V.[1] + dti_bj_12 t))
        (fun t -> fun res -> 
            res.[0] <- dti_bj_00 t
            res.[1] <- dti_bj_11 t)
        50 0 100 (Array.zeroCreate 2)