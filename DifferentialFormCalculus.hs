module DifferentialFormCalculus where

import Data.List

data Sign = Plus | Minus | Zero deriving (Show, Eq)

findPermutationSign :: (Ord a) => [a] -> Int
findPermutationSign = findPermutationSignBy compare

findPermutationSignBy :: (Eq a) => (a -> a -> Ordering) -> [a] -> Int
findPermutationSignBy ord x = if res == -1 then 0
  else if even res then 1
    else (-1)
  where res = countTranspositionsBy ord x

countTranspositionsBy :: (Eq a) => (a -> a -> Ordering) -> [a] -> Int -- terribly written function :(
countTranspositionsBy _ [] = 0
countTranspositionsBy f q@(x:xs) = if areThereDups then -1
  else if smallest == x then countTranspositionsBy f xs
    else ((countTranspositionsBy f (swap x q)) + 1) where
      areThereDups = (length q) /= (length qSorted)
      qSorted = nub $ sortBy f q
      smallest = head qSorted
      swap x xs = beginning++rest where
        (l1,l2) = splitAt pos xs 
        beginning = (head l2):(tail l1)
        rest = x:(tail l2)
        pos = head (findIndices (x ==) qSorted)


uniqStrIncNtuple :: Int -> [a] -> [[a]]
uniqStrIncNtuple 0 _      = [[]]
uniqStrIncNtuple _ []     = []
uniqStrIncNtuple n (x:xs) = (map (x:) $ uniqStrIncNtuple (n-1) xs) ++ uniqStrIncNtuple n xs

ngenFromBasis :: Int -> [a] -> [[a]]
ngenFromBasis n b = uniqStrIncNtuple n b

formComps :: (Num a, Eq a) => Int -> [bvs] -> ([bvs] -> a) -> [([bvs], a)]
formComps n bvs f = filter ((/=0) . snd ) . zip (ngenFromBasis n bvs) $ map f (ngenFromBasis n bvs)

alternate :: (Eq vs, Fractional cn) => ([vs] -> cn) -> [vs] -> cn
alternate f l = (divBy *) . sum . zipWith (\a b -> fromIntegral a * b) sgns $ map f perms
  where
    divBy = recip . fromIntegral $ product [1..(length l)]
    perms = permutations l
    sgns = map (findPermutationSignBy cmp) perms
    cmp a b = compare (lkupFn a) (lkupFn b)
    lkupFn = flip lookup (zip l [0..])

intPr :: v -> ([v] -> a) -> [v] -> a
intPr x formOnList = formOnList . (x:)

wdgPrOnComps :: (Eq a, Num a) => [(String, a)] -> [(String, a)] -> [(String, a)]
wdgPrOnComps f1 f2 = collectLikeTerms $ nonZeroes [combine c1 c2 | c1 <- f1, c2 <- f2] where
  combine (v1,n1) (v2,n2) = (sort v12, fromIntegral s * n1*n2) where
    s = findPermutationSign v12
    v12 = v1 ++ v2
  nonZeroes = filter ((/=0) . snd)
  collectLikeTerms = map reduceGroups . groupBy (\a b-> fst a == fst b) . sortOn fst
  reduceGroups l = (head ns, sum vals) where
    (ns, vals) = unzip l
