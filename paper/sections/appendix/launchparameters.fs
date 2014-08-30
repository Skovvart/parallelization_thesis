let launchParameters iterations = 
    if iterations <= 32 then 1, iterations
    else if iterations <= 14*32 then divup iterations 32, 32
    else
        let rec optParams blocks threads =
            if blocks * threads >= iterations then blocks, threads
            else
                let dBlocks = blocks + 28 // optimal solution not guaranteed
                let dThreads = min 512 (threads + 32)
                if blocks * dThreads >= iterations then blocks, dThreads
                else if dBlocks * dThreads >= iterations then 
                    let dBlocks = [blocks .. dBlocks] |> List.find (fun b -> b * dThreads >= iterations)
                    dBlocks, dThreads
                else
                    if dThreads = 512 then optParams dBlocks 512
                    else optParams blocks dThreads
        optParams 14 32