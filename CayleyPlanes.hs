module CayleyPlanes where

import ComplexOctonions
import ComplexNumbers (ComplexNumber, i)
import Data.List
import DifferentialFormCalculus

gen3VsIn7d :: [[Octonion]]
gen3VsIn8d :: [[Octonion]]
gen4VsIn7d :: [[Octonion]]
gen4VsIn8d :: [[Octonion]]
gen3VsIn7d = ngenFromBasis 3 imagOctGens
gen3VsIn8d = ngenFromBasis 3 octGens
gen4VsIn7d = ngenFromBasis 4 imagOctGens
gen4VsIn8d = ngenFromBasis 4 octGens

vName :: Octonion -> Char
vName = unjust . flip lookup (zip octGens "01234567")
  where
    unjust (Just v) = v
    unjust Nothing  = error "non generator vector being named!"

vNames :: [Octonion] -> String
vNames = map vName

nvecNames :: [Octonion] -> String
ncovNames :: [Octonion] -> String
nvecNames os = "e_{" ++ vNames os ++ "}"
ncovNames os = "e^{" ++ vNames os ++ "}"

removeVecCov = init . drop 3

formComps7 n f = map (\(a,b) -> (vNames a, b)) $ formComps n imagOctGens f
formComps8 n f = map (\(a,b) -> (vNames a, b)) $ formComps n octGens f
formComps7' n f = map (\(a,b) -> (ncovNames a, b)) $ formComps n imagOctGens f
formComps8' n f = map (\(a,b) -> (ncovNames a, b)) $ formComps n octGens f

crProd :: Octonion -> Octonion -> Octonion
crProd u v = imaginary ((conj v) * u)

crProdOnList :: [Octonion] -> Octonion
crProdOnList [u,v] = crProd u v

phi :: Octonion -> Octonion -> Octonion -> ComplexNumber
phi u v w = dotProd (crProd u v) w

phiOnList :: [Octonion] -> ComplexNumber
phiOnList [u,v,w] = phi u v w
phiOnList _ = error "phi requires exactly 3 vectors"

phiComps = formComps7 3 phiOnList

contractPhi x = formComps7 2 (intPr x phiOnList)

metricAndVol x = wdgPrOnComps (contractPhi x) $ wdgPrOnComps (contractPhi x) phiComps

-- we want to express triCrProd as alternation
-- of u(\bar v w).

preTrCrPrOnList = alternate f where
  f [u, v, w] = (u * (conj v)) * w

-- this definition of triple cross product agrees with the one below.

triCrProd :: Octonion -> Octonion -> Octonion -> Octonion
triCrProd u v w = (scalarMult 0.5) ((u * cnjv) * w - (w * cnjv) * u) where cnjv = conj v

triCrProdOnList :: [Octonion] -> Octonion
triCrProdOnList [u,v,w] = triCrProd u v w
triCrProdOnList _ = error "triCrProd requires exactly 3 vectors"

capPhi :: Octonion -> Octonion -> Octonion -> Octonion -> ComplexNumber
capPhi x u v w = dotProd x (triCrProd u v w)

capPhiOnList :: [Octonion] -> ComplexNumber
capPhiOnList [x,u,v,w] = capPhi x u v w
capPhiOnList _ = error "capPhi is a 4 form! It requires 4 vectors to be evaluated."

capPhiComps = formComps8 4 capPhiOnList

latexPrintCvalForms :: [(String, ComplexNumber)] -> String
latexPrintCvalForms = concatMap putCoefficient
  where
    putCoefficient (s,1) = '+':s
    putCoefficient (s,-1) = '-':s
    putCoefficient (s,x) = (show x)++s

-- this function is not convenient yet.
latexPrintVvalForms :: [(String, Octonion)] -> String
latexPrintVvalForms = concatMap putCoefficient
  where
    putCoefficient (s,x) = " + " ++ (show x) ++ s

latexPrintCapPhi = latexPrintCvalForms $ formComps8' 4 capPhiOnList

preQuadCrProdOnList = alternate f
  where
    f [x,u,v,w] = negate $ ((x * conj u) * v) * conj w
    --f [x,u,v,w] = conj x*triCrProdOnList [u,v,w]

quadCrPr x u v w = scalarMult (-0.25) (
  (triCrProd x u v) * conj w
  -(triCrProd w x u) * conj v
  +(triCrProd v w x) * conj u
  -(triCrProd u v w) * conj x)

quadCrPrOnList [x,u,v,w] = quadCrPr x u v w

isAllOnes = map (\q -> ((capPhiOnList q)**2) + normSquared (imaginary $ preQuadCrProdOnList q)) gen4VsIn8d

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
