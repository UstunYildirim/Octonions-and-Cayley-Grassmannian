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
    
