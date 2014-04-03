open System.IO

let getValues (line:string) =
    let line = (line.Substring 4).Trim()
    let values = line.Split ' ' |> Array.filter (fun v -> v <> "")
    Seq.map (fun (s:string) -> s.Replace (',', '.') |> float) values

let compareValues (fst:seq<float>) snd tolerance =
    if Seq.length fst <> Seq.length snd then
        printfn "Error: The two sequences do not have the same length"
        Seq.empty 
    else
        Seq.map2 (fun f -> fun s -> if abs (f-s) > tolerance then false else true) fst snd

let countErrors results = 
    Seq.fold 
        (fun errCount -> fun rowRes -> errCount + (Seq.fold (fun rowErrCount -> fun res -> rowErrCount + if res then 0 else 1)
        0 rowRes)) 0 results

let skipEmptyLines lines =
    Array.filter (fun line -> line <> "") lines


//results/double-C-sharp.txt results/double-cuBase.txt
//results/double-C-sharp.txt results/float-C-Sharp.txt
//results/double-cuBase.txt results/float-cuBase.txt
//results/float-C-sharp.txt results/float-cuBase.txt
//results/float-C-sharp.txt results/float-cuda.txt
//results/float-cuBase.txt results/float-cuda.txt


let compare fileNames errorTolerance =
    let filesToCompare = Array.map (fun fileName -> File.ReadLines fileName |> Seq.toArray) fileNames
    let f1 = skipEmptyLines filesToCompare.[0] |> Array.map getValues
    let f2 = skipEmptyLines filesToCompare.[1] |> Array.map getValues
    let testResults = Array.map2 (fun f -> fun s -> compareValues f s errorTolerance) f1 f2
    let resultLength = Seq.fold (fun count -> fun (e:seq<bool>) -> count + (Seq.length e)) 0 testResults
    let totalErrors = countErrors testResults
    printfn "%d values compared, %d total errors with error tolerance %.0e" resultLength totalErrors errorTolerance

let compareRange (fileNames:string[]) = 
    printfn "Comparing %s with %s" fileNames.[0] fileNames.[1]
    compare fileNames 0.0
    compare fileNames 1e-15
    compare fileNames 1e-6
    compare fileNames 1e-5
    compare fileNames 1e-4
    printfn ""

[<EntryPoint>]
let main argv = 
    compareRange [|"double-C-sharp.txt"; "double-cuBase.txt"|]
    compareRange [|"double-C-sharp.txt"; "float-C-sharp.txt"|]
    compareRange [|"double-cuBase.txt"; "float-cuBase.txt"|]
    compareRange [|"float-C-sharp.txt"; "float-cuBase.txt"|]
    compareRange [|"float-C-sharp.txt"; "float-cuda.txt"|]
    compareRange [|"float-cuBase.txt"; "float-cuda.txt"|]

    0 // return an integer exit code
