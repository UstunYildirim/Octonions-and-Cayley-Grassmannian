module Jacobians where

import Data.List
import Data.Maybe
import ComplexNumbers (ComplexNumber, i)
import ComplexOctonions
import LaurentPolynomials
import DifferentialFormCalculus
import MatrixCalculus
import CayleyPlanes
import MaximalTorus
import EigenvalueList
import Plucker

normTransDefEqsAsPoly = map (Poly . map createTerm) normTransDefEqs
  where
    createTerm (var,coeff) = Term coeff [(var,1)]

localizeNewDefEqs locChart = map (simplify . localizeTo locChart) normTransDefEqsAsPoly

latexNewLocalizedDefEqs locChart = map (mapVars typeAsNewPlucker) $ localizeNewDefEqs locChart
  where
    typeAsNewPlucker = ("\\wil q_{"++) . (++"}")

-- we will compute the Jacobians at the origin
-- this means the quadratic terms will not contribute
-- but Plucker relations give at least quadratic
-- replacements. so unless a variable is already
-- in a local variable it will not contribute.
-- we use this shortcut to compute Jacobians
localVars locChart = tail . take 17 $ sortOn (distance locChart) allCoordinates

jacobianOfDefEqsAtO :: String -> Matrix (LPoly ComplexNumber)
jacobianOfDefEqsAtO locChart = Matrix $ map createRow normTransDefEqs
  where
    createRow ls = map lookupEntry (localVars locChart) where
      lookupEntry lvar = let val = lookup lvar ls in
        if isNothing val then 0 else fromJust val

genKerOfJac locChart = map snd zeroColsWvecs ++ generatorsFromNonZeroPart
  where
    allCols = transpose . matrixAsList $ jacobianOfDefEqsAtO locChart
    allInfo = zip allCols (localVars locChart)
    (zeroColsWvecs, nonZeroColsWvecs) = partition test allInfo
    test = (== zeroRow) . fst
    zeroRow = head (matrixAsList (zeroMatrix 1 7))
    nonZeroCols = map fst nonZeroColsWvecs
    generatorsFromNonZeroPart = map createVectors $ matchCols nonZeroCols
    createVectors (c1,c2) = if c1 == c2
      then fromJust (lookup c1 allInfo) ++ " - " ++ fromJust (lookup c2 allInfo)
      else fromJust (lookup c1 allInfo) ++ " + " ++ fromJust (lookup c2 allInfo)
    matchCols [] = []
    matchCols (c:cs) = if isNothing mayCol
      then matchCols cs
      else 
        if eigValLookUp (fromJust $ lookup c allInfo) == eigValLookUp (fromJust $ lookup col allInfo)
          then (c,col):matchCols (delete col cs) 
          else error $ "Eigenvalues do not match! " ++
            "(Actually this error might be eliminated but " ++
            "implementation is not necessary at the moment.)"
      where
        mayCol = find (doColsPair c) cs
        col = fromJust mayCol
    doColsPair c1 c2 = isColsLinDep c1 c2 && doColsHaveSameEigVal c1 c2
    isColsLinDep c1 c2 = c1 == c2 || (map (*(-1)) c1) == c2
    doColsHaveSameEigVal c1 c2 =
      eigValLookUp (fromJust $ lookup c1 allInfo) == eigValLookUp (fromJust $ lookup c2 allInfo)

