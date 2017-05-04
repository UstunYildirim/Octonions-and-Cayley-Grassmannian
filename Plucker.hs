module Plucker where

import Data.List
import Data.Maybe
import ComplexOctonions
import DifferentialFormCalculus
import CayleyPlanes
import LaurentPolynomials


-- example of a plucker relation: 
-- p_{xyza}p_{bcde}
--  = p_{xyzb}p_{acde} - p_{xyzc}p_{abde} + p_{xyzd}p_{abce} - p_{xyze}p_{abcd}
-- the function below takes "xyza" "bcde" as input
-- and creates a LPoly with variables of the same form
-- note: this works for arbitrary length (not only 4 chars)
genPluckerRel :: (Num a) => String -> String -> LPoly a
genPluckerRel xyza bcde = res where
  xyz = init xyza
  abcde = last (xyza) : bcde
  ln = length abcde
  diffCuts = map (flip splitAt abcde) [1..ln-1]
  variations = map (\(ab,cde)->(xyz ++ [head cde], ab ++ tail cde)) diffCuts
  res = Poly [createTerm package | package <- zip variations (cycle [1,-1])]
  createTerm ((v1,v2),c) = Term c [(v1,1), (v2,1)]

-- same as above except this one 
-- reorders variables and set
-- coefficient equal to 0 if there
-- are repeated indices
genOrderedPluckerRel :: (Num a) => String -> String -> LPoly a
genOrderedPluckerRel xyza bcde = res where
  xyz = init xyza
  abcde = last (xyza) : bcde
  ln = length abcde
  diffCuts = map (flip splitAt abcde) [1..ln-1]
  variations = map (\(ab,cde)->(xyz ++ [head cde], ab ++ tail cde)) diffCuts
  res = Poly [createOrderedTerm package | package <- zip variations (cycle [1,-1])]
  createOrderedTerm ((v1,v2),c) = Term c' [(v1',1), (v2',1)] where
    c1 = findPermutationSign v1
    c2 = findPermutationSign v2
    c' = c*c1*c2
    v1' = sort v1
    v2' = sort v2

genOrderedNonTrivPluckerRel :: (Eq a, Num a) => String -> String -> LPoly a
genOrderedNonTrivPluckerRel xyza abcd = correctCoeff $ genOrderedPluckerRel chosenPerm abcd where
  correctCoeff = if even permCoeff then id else ((Poly [Term (fromInteger (-1)) []]) *)
  cyclicPerms = [take len $ drop i $ cycle xyza | i <- [0..len-1] ]
  chosenPerm = head $ ((filter isNonTriv cyclicPerms) ++ [xyza])
  permCoeff = fromJust $ findIndex (chosenPerm ==) cyclicPerms
  isNonTriv perm = not $ elem (last perm) abcd
  len = length xyza

allCoordinates = ngenFromBasis 4 "01234567"

distance :: String -> String -> Int
distance c1 c2 = 4 - length (intersect c1 c2)

localizeTo :: (Eq a, Num a) => String -> LPoly a -> LPoly a
localizeTo locChart = reduceToLocalVars where
  allNonLocalVarsSorted = drop 17 $ sortOn (distance locChart) allCoordinates
  replaceVar var = plugInIdentity var (genOrderedNonTrivPluckerRel var locChart)
  reduceToLocalVars = foldr (.) id $ map replaceVar allNonLocalVarsSorted

localizeDefEqs locChart = map (simplify . localizeTo locChart) definingEquations

latexLocalizedDefEqs locChart = map (mapVars typeAsPlucker) $ localizeDefEqs locChart
  where
    typeAsPlucker = ("p_{"++) . (++"}")


