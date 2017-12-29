module Rand

/// A probability (subset of unit interval [0,1]).
type Prob = float

/// A probability distribution for sample space 'a.
type Dist<'a> = 
    | D of seq<'a * Prob>

[<AutoOpen>]
module private Util =

    let cartesian (a:seq<'a>) (b:seq<'b>) : seq<'a * 'b> =
      Seq.collect (fun a -> b |> Seq.map (fun b -> a,b)) a

    let rec factorial (n:int) : int =
        if n <= 0 then 1
        else n * factorial (n - 1)

    let binomialCoef (n:int) (k:int) : int = 
        (factorial n) / (factorial (n - k) * factorial k)

/// Operations on probability distributions.
module Dist =

    let private unD (D s) = s

    let map (f:'a -> 'b) (d:Dist<'a>) : Dist<'b> =
        unD d |> Seq.map (fun (a,p) -> f a,p) |> D

    let scale (xs:seq<'a * Prob>) : Dist<'a> =
        let q = xs |> Seq.sumBy snd in
        D (xs |> Seq.map (fun (a,p) -> a, p/q))

    let zipWith (f:'a -> 'b -> 'c) (da:Dist<'a>) (db:Dist<'b>) : Dist<'c> =
        cartesian (unD da) (unD db) 
        |> Seq.map (fun ((a,p),(b,q)) -> f a b, p * q)
        |> D
    
    let joint (da:Dist<'a>) (db:Dist<'b>) : Dist<'a * 'b> =
        zipWith (fun a b -> a,b) da db

    let certainly (a:'a) : Dist<'a> =
        D <| Seq.singleton (a,1.0)

    let impossible<'a> : Dist<'a> = 
        D Seq.empty

    let bind (f:'a -> Dist<'b>) (da:Dist<'a>) : Dist<'b> =
        D <| seq { 
          for (a,p) in (unD da) do 
            for (b,q) in (unD (f a)) -> (b, p * q) }

    let private fromPmfInt (pmf:int -> Prob) : Dist<int> =
        D <| Seq.initInfinite (fun k -> k, pmf k)

    /// The Binomial distribution models the number of successes in N trials.
    /// The support of the distribution is the number of successes in N trials {0...n}.
    let binomial (n:int) (p:Prob) : Dist<int> =
        let pmf k = float (binomialCoef n k) * (pown p k) * (pown (1.0 - p) (n - k))
        fromPmfInt pmf

    /// The geometric distribution models the number of Bernoulli trials required for success.
    /// The support is a number of trials {0,....}.
    /// n - the number of trials
    /// p - probability of success
    let geometric (p:Prob) : Dist<int> =
        let pmf k = (pown (1.0 - p) (k - 1)) * p
        fromPmfInt pmf

    /// The Poisson distribution models the number of occurrences of events during a time interval.
    /// Events occur independently, with a constant rate.
    /// rate - the average rate
    let poisson (rate:float) : Dist<int> =
        let pmf k = (pown rate k) * (exp -rate) / float (factorial k)
        fromPmfInt pmf
    
    let uniform (xs:'a list) : Dist<'a> =
        let p = 1.0 / float (List.length xs)
        xs |> Seq.map (fun a -> a,p) |> D

    