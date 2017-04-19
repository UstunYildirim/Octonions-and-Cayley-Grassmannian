module CayleyPlanes where

import ComplexOctonions
import ComplexNumbers (ComplexNumber, i)
import Data.List

uniqStrIncNtuple :: Int -> [a] -> [[a]]
uniqStrIncNtuple 0 _      = [[]]
uniqStrIncNtuple _ []     = []
uniqStrIncNtuple n (x:xs) = (map (x:) $ uniqStrIncNtuple (n-1) xs) ++ uniqStrIncNtuple n xs

gen3VsIn7d :: [[Octonion]]
gen3VsIn8d :: [[Octonion]]
gen4VsIn7d :: [[Octonion]]
gen4VsIn8d :: [[Octonion]]
gen3VsIn7d = uniqStrIncNtuple 3 imagOctGens
gen3VsIn8d = uniqStrIncNtuple 3 octGens
gen4VsIn7d = uniqStrIncNtuple 4 imagOctGens
gen4VsIn8d = uniqStrIncNtuple 4 octGens

nNames37 = uniqStrIncNtuple 3 "1234567"
nNames38 = uniqStrIncNtuple 3 "01234567"
nNames47 = uniqStrIncNtuple 4 "1234567"
nNames48 = uniqStrIncNtuple 4 "01234567"

nvecNames37 = map (("e_{" ++) . (++ "}")) nNames37
nvecNames38 = map (("e_{" ++) . (++ "}")) nNames38
nvecNames47 = map (("e_{" ++) . (++ "}")) nNames47
nvecNames48 = map (("e_{" ++) . (++ "}")) nNames48

nCovNames37 = map (("e^{" ++) . (++ "}")) nNames37
nCovNames38 = map (("e^{" ++) . (++ "}")) nNames38
nCovNames47 = map (("e^{" ++) . (++ "}")) nNames47
nCovNames48 = map (("e^{" ++) . (++ "}")) nNames48

formComps37 f = filter ((/=0) . snd ) . zip nCovNames37 $ map f gen3VsIn7d
formComps38 f = filter ((/=0) . snd ) . zip nCovNames38 $ map f gen3VsIn8d
formComps47 f = filter ((/=0) . snd ) . zip nCovNames47 $ map f gen4VsIn7d
formComps48 f = filter ((/=0) . snd ) . zip nCovNames48 $ map f gen4VsIn8d

crProd :: Octonion -> Octonion -> Octonion
crProd x y = imaginary (x * (conj y))

phi :: Octonion -> Octonion -> Octonion -> ComplexNumber
phi x y z = dotProd x (crProd y z)

phiOnList :: [Octonion] -> ComplexNumber
phiOnList [x,y,z] = phi x y z
phiOnList _ = error "phi requires exactly 3 vectors"

-- we want to express triCrProd as alternation
-- of u(\bar v w).

preTrCrPr u v w = scalarMult (1/6) $ f u v w - f u w v
                                    + f w u v - f w v u
                                    + f v w u - f v u w
  where
    f u v w = u * ((conj v) * w)

preTrCrPrOnList [u,v,w] = preTrCrPr u v w

-- this definition of triple cross product agrees with the one below.

triCrProd :: Octonion -> Octonion -> Octonion -> Octonion
triCrProd x y z = (scalarMult 0.5) (x * (cnjy * z) - z * (cnjy * x)) where cnjy = conj y

triCrProdOnList :: [Octonion] -> Octonion
triCrProdOnList [x,y,z] = triCrProd x y z
triCrProdOnList _ = error "triCrProd requires exactly 3 vectors"

capPhi :: Octonion -> Octonion -> Octonion -> Octonion -> ComplexNumber
capPhi x y z w = dotProd x (triCrProd y z w)

capPhiOnList :: [Octonion] -> ComplexNumber
capPhiOnList [x,y,z,w] = capPhi x y z w
capPhiOnList _ = error "capPhi is a 4 form! It requires 4 vectors to be evaluated."

latexPrintCapPhi = concatMap putPM $ formComps48 capPhiOnList
  where
    putPM (s,1) = '+':s
    putPM (s,-1) = '-':s

-- 
-- commutator :: Octonion -> Octonion -> Octonion
-- commutator x y = x*y - y*x
-- associator :: Octonion -> Octonion -> Octonion -> Octonion
-- associator x y z = (x * y) * z - x * (y * z)
-- associatorOnList :: [Octonion] -> Octonion
-- associatorOnList [x,y,z] = associator x y z
-- associatorOnList _ = error "associator requires exactly 3 vectors"
--
-- coassociator :: Octonion -> Octonion -> Octonion -> Octonion -> Octonion
-- coassociator x y z w = (-2)*((sm (dp y' (z'*w')) x')+(sm (dp z' (x'*w')) y')+(sm (dp x' (y'*w')) z')+(sm (dp y' (x'*z')) w')) where
--                          [x',y',z',w'] = map imaginary [x,y,z,w]
--                          sm = scalarMult
--                          dp = dotProd
-- -- End of NOTE1
-- 
-- 
-- -- psi is a 4-form on R^7
-- -- it is equal to starPhi so the definition below is just for reference
-- psi :: Octonion -> Octonion -> Octonion -> Octonion -> ComplexNumber 
-- psi x y z w = (1/2) * (dotProd x (associator y z w))
-- psiOnList :: [Octonion] -> ComplexNumber
-- psiOnList [x,y,z,w] = psi x y z w
-- psiOnList _ = error "psi requires exactly 4 vectors."
-- 
-- -- dual form of phi; starPhi is a 4-form on R^7
-- starPhi :: Octonion -> Octonion -> Octonion -> Octonion -> ComplexNumber
-- starPhi x y z w = dotProd x (tripleCrProd y z w)
-- starPhiOnList :: [Octonion] -> ComplexNumber
-- starPhiOnList [x,y,z,w] = starPhi x y z w
-- starPhiOnList _ = error "starPhi requires exactly 4 vectors"
-- 
-- 
-- -- chi is a vector (imaginary octonion) valued 3-form on R^7
-- -- it is obtained by "lowering one of the indices of *phi"
-- -- i.e. *phi(x,y,z,w) = <chi(x,y,z),w>
-- -- chi may also be computed as imaginary.tripleCrProd
-- -- or (/3).jacob' where jacob' [x,y,z] is the jacobi identity computed on x,y,z with cross prod
-- chi :: Octonion -> Octonion -> Octonion -> Octonion
-- chi x y z = foldr (\u v->(f u)+v) 0 imagOctGens where
--   f o = scalarMult ((starPhi x y z o)/(normSquared o)) o 
-- chiOnList :: [Octonion] -> Octonion
-- chiOnList [x,y,z] = chi x y z
-- chiOnList _ = error "chi requires exactly 3 vectors"
-- -- End of NOTE2
-- 
-- 
-- -- auxiliary 
-- twiceIm4foldCrProd :: Octonion -> Octonion -> Octonion -> Octonion -> Octonion
-- twiceIm4foldCrProd x y z w = (coassociator x y z w) + (sm x1 (associator y z w)) + (sm y1 (associator z x w)) + (sm z1 (associator x y w)) + (sm w1 (associator y x z)) where
--                                   [x1,y1,z1,w1] = map real [x,y,z,w]
--                                   sm = scalarMult
-- 
-- -- As its name suggests, imaginary part of 4 fold cross product is R^7-valued 4 form on R^8
-- im4foldCrProd :: Octonion -> Octonion -> Octonion -> Octonion -> Octonion
-- im4foldCrProd x y z w = scalarMult 0.5 $ twiceIm4foldCrProd x y z w
-- im4foldCrProdOnList :: [Octonion] -> Octonion
-- im4foldCrProdOnList [x,y,z,w] = im4foldCrProd x y z w
-- im4foldCrProdOnList _ = error "im4foldCrProdOnList requires exactly 4 vectors"
-- 
-- 
-- pullback :: (a -> b) -> ([b] -> c) -> [a] -> c
-- pullback function form input = form (map function input)
-- 
-- pairNonZeroResults                :: (Num a, Eq a) => ([Octonion] -> a) -> [[Octonion]] -> [(String,a)]
-- pairNonZeroResultsForIm4FldCrProd :: [(String, Octonion)]
-- pairNonZeroResultsForChi          :: [(String, Octonion)]
-- pairNonZeroResultsForPhi          :: [(String, ComplexNumber)]
-- pairNonZeroResultsForStarPhi      :: [(String, ComplexNumber)]
-- pairNonZeroResultsForCapPhi       :: [(String, ComplexNumber)]
-- 
-- pairNonZeroResults form listOfGens  = filter (\(_,y)->(y/=0)) $ zip (map printShortNames listOfGens) $ map form listOfGens
-- pairNonZeroResultsForIm4FldCrProd   = pairNonZeroResults im4foldCrProdOnList genOf4VectsIn8dim
-- pairNonZeroResultsForChi            = pairNonZeroResults chiOnList genOf3VectsIn7dim
-- pairNonZeroResultsForPhi            = pairNonZeroResults phiOnList genOf3VectsIn7dim
-- pairNonZeroResultsForStarPhi        = pairNonZeroResults starPhiOnList genOf4VectsIn7dim
-- pairNonZeroResultsForCapPhi         = pairNonZeroResults capPhiOnList genOf4VectsIn8dim
-- 
-- printShortNames :: [Octonion] -> String -- prints of the form e^{148}
-- printShortNames = (++"}").("e^{"++).foldr f "" where
--   f e st = if (r == Nothing) then st else (s r) ++ st where
--     r          = Data.List.elemIndex e octGens
--     s (Just a) = show a
--     s Nothing  = ""
-- 
-- printShortName :: Octonion -> String -- prints of the form e_4, -e_3
-- printShortName x = text where
--     text = if (n==Nothing) then if (n'==Nothing) then "" else '-':(get n') else get n where
--      get Nothing  = ""
--      get (Just a) = "e_"++(show a)
--      n            = Data.List.elemIndex x octGens
--      n'           = Data.List.elemIndex (-x) octGens
-- 
-- evaluateFormOnGens :: (Ord x, Num x) => ([Octonion]->x) -> [[Octonion]] -> [([String],x)]
-- evaluateFormOnGens form gens = map refineList groupedResults where
--   refineList l   = ((map fst l),snd (l!!0))
--   groupedResults = groupBy (\(_,x) (_,y)-> x == y) sortedResults where
--     sortedResults = sortOn snd nonZeroPairedResults where
--       nonZeroPairedResults = filter (\(_,x)->x/=0) $ zip prettyGens results where
--         prettyGens = map printShortNames gens
--         results    = map form gens
-- 
-- pickCoefficients :: Eq a => [(a1,a)] -> a -> [(a1,a)]
-- pickCoefficients pairedResults coeff = filter (\x->(snd x)==coeff) pairedResults
-- 
-- prettyPrintForm :: [(String,Octonion)] -> [Octonion] -> String
-- prettyPrintForm pairedResults coeffList = intercalate " \\\\\n+" $ map (prettyPrintOneComponent pairedResults) coeffList where
--   prettyPrintOneComponent pairedResults' coeff = "("++(concatMap ("-"++) . map fst $ pickCoefficients pairedResults' (-coeff))++(concatMap ("+"++) . map fst $ pickCoefficients pairedResults' coeff)++")"++printShortName coeff
-- 
-- prettyPrintIm4fldCrProd :: String
-- prettyPrintIm4fldCrProd = prettyPrintForm pairNonZeroResultsForIm4FldCrProd imagOctGens --Latex Ready
-- 
-- printLineByLine :: (Show a) => [a] -> IO()
-- printLineByLine = mapM_ (putStrLn.show)
-- 
-- 
-- -- TEST CODE (for associative planes)
-- --
-- 
-- -- these can be refactored thru uniqStrIncNtuple
-- uniqStrIncPairs :: [a] -> [(a,a)]
-- uniqStrIncPairs l = [(x,y) | (x:ys) <- tails l, y <- ys]
-- uniqStrIncTriples :: [a] -> [(a,a,a)]
-- uniqStrIncTriples l = [(x,y,z) | (x:ys) <- tails l, (y:zs) <- tails ys, z<-zs]
-- uniqStrIncQuadruples :: [a] -> [(a,a,a,a)]
-- uniqStrIncQuadruples l = [(x,y,z,w) | (x:ys)<-tails l, (y:zs)<-tails ys, (z:ws)<-tails zs, w<-ws]
-- -- END OF these can be refactored thru uniqStrIncNtuple
-- 
-- someComplexOctonions = nub . sort $ filter (/=0) [a*b+c*d | a<-[1,ii], c<-[1,ii,-1,-ii],(b,d)<-uniqStrIncPairs octGens]
-- somePureImagComplexOctonions = nub . sort $ filter (/=0) [a*b+c*d | a<-[1,ii,-1,ii], c<-[1,ii,-1,-ii], (b,d)<-uniqStrIncPairs imagOctGens]
-- 
-- (zeroes,nonZeroes) = partition ((==0).normSquared) someComplexOctonions
-- (pureImagZeroes,pureImagNonZeroes) = partition ((==0).normSquared) somePureImagComplexOctonions
-- 
-- units = map makeUnit nonZeroes
-- pureImagUnits = map makeUnit pureImagNonZeroes
-- 
-- isLinIndep :: (Octonion,Octonion)->Bool
-- isLinIndep = not . isReal . uncurry (/)
-- 
-- nonTrivPairs = filter isLinIndep . uniqStrIncPairs -- assumes distinct elements
-- 
-- linIndepNonZeros = nonTrivPairs nonZeroes 
-- pureImagLinIndepNonZeros = nonTrivPairs pureImagNonZeroes 
-- linIndepUnits = nonTrivPairs units
-- pureImagLinIndepUnits = nonTrivPairs pureImagUnits
-- 
-- givesAssocPlane :: (Octonion, Octonion)->Bool -- assumes unit vectors (norm = 1)
-- givesAssocPlane (a,b) = let v = crProd a b in normSquared v == phi a b v
-- 
-- isProbPlane :: (Octonion, Octonion)->Bool -- assumes 0 normed vectors
-- isProbPlane (a,b) = let v = crProd a b in if normSquared v /= 0 then False else 0/=phi a b v
-- 
-- rebelPlanes = filter (not . givesAssocPlane) pureImagLinIndepUnits -- there are no rebels
-- zeroRebels = filter isProbPlane $ nonTrivPairs pureImagZeroes -- there are no rebels
-- 
-- v1 = ii*oi+oj
-- v2 = ok-ol
-- v3 = crProd v1 v2
-- 
-- cmp' l = (phiOnList l)**2 == 1
-- reducedImZeroes = nub . sort $ filter ((==0).normSquared) $ filter (/=0) [b+c*d | c<-[1,ii,-1,-ii],(b,d)<-uniqStrIncPairs imagOctGens]
-- lotsOf3Gens = uniqStrIncNtuple 3 reducedImZeroes
-- pff' = head $ filter isOrthogonalList $ filter cmp' lotsOf3Gens
-- 
-- -- END OF TEST CODE (for associative planes)
-- --
-- -- TEST CODE (for Cayley Planes)
-- 
-- reducedZeroes = nub . sort $ filter ((==0).normSquared) $ filter (/=0) [b+c*d | c<-[1,ii,-1,-ii],(b,d)<-uniqStrIncPairs octGens]
-- lotsOf4Gens = uniqStrIncNtuple 4 reducedZeroes
-- 
-- defn = map ((**2).capPhiOnList) lotsOf4Gens
-- lemma = map (normSquared.im4foldCrProdOnList) lotsOf4Gens
-- 
-- cmp l = (capPhiOnList l)**2 + normSquared (im4foldCrProdOnList l) == 1
-- 
-- vs = head $ filter cmp lotsOf4Gens
-- 
-- pairwiseDistinctBy :: [a] -> (a->a->Bool) -> [a] -- assumes the comparison function is commutative
-- pairwiseDistinctBy (x:xs) isDistinct = if distinct then (x:rest) else rest where
--   distinct = and $ map (isDistinct x) xs
--   rest = pairwiseDistinctBy xs isDistinct
-- pairwiseDistinctBy _ _ = []
-- 
-- isOrthogonalList :: [Octonion] -> Bool
-- isOrthogonalList [] = True
-- isOrthogonalList (x:xs) = (and $ (map ((==0).(dotProd x)) xs)) && (isOrthogonalList xs)
-- 
-- pff = head $ filter isOrthogonalList $ filter cmp lotsOf4Gens
-- 
-- bunchOfImagOcts = uniqStrIncNtuple 3 imagOctGens
-- jacob op [x,y,z] = (op x (op y z)) + (op y (op z x)) + (op z (op x y))
-- jacob' = (/3).jacob crProd -- R^7 valued 3 form on R^7
-- 
-- pairNonZeroResultsForJacob' = pairNonZeroResults jacob' genOf3VectsIn7dim
-- 
-- pbsByRightMult = map (\x->pullback (*x) capPhiOnList) (octGens ++ map (*ii) octGens)
-- theList = nub $ map (\x->pairNonZeroResults x genOf4VectsIn8dim) pbsByRightMult
-- theList' = pairwiseDistinctBy theList (\x y-> not . and $ zipWith sndCorClose x y)
-- 
-- sndCorClose = (\(_,x) (_,y)-> abs (x-y) < 1e-10)
-- 
-- -- END OF TEST CODE (for Cayley Planes)
-- --}
