module Test where

import Data.List
import ComplexNumbers (ComplexNumber, i)
import ComplexOctonions
import CayleyPlanes
import DifferentialFormCalculus
import MatrixCalculus
import LaurentPolynomials
import MaximalTorus
import EigenvalueList
import Plucker
import Jacobians

-- nulls = [1+ii*oi, oi+ii*oj, oj + ii*ok]
-- 
-- isNormEqTrueForQuadCrPr = map (normSquared . quadCrPrOnList) gen4VsIn8d == replicate (length gen4VsIn8d) 1
-- 
-- capPhiSquared = wdgPrOnComps capPhiComps capPhiComps
-- 
-- sumOfCapPhiAndImQd = map (\q -> ((capPhiOnList q)**2) + normSquared (imaginary $ preQuadCrProdOnList q)) gen4VsIn8d
-- 
-- isAllOnes :: [ComplexNumber] -> Bool
-- isAllOnes [] = True
-- isAllOnes (1:xs) = isAllOnes xs
-- isAllOnes _ = False

cxToPoly c = Poly [Term c []]

imQdCrPrAsMat' = transposeMat . fmap cxToPoly $ Matrix [octonionAsCxList $ imQuadCrPrOnList l | l <- basis ] where
  basis = uniqStrIncNtuple 4 octGens

-- below is an attempt to compute derivative at \wil "0246" 
-- without referring to transformed equations:
-- what we do is take wedge of the four vectors 
-- i.e. flip indMatOnExtAlg 4 of the 8x4 matrix
-- whose columns are these vectors in the standard basis
-- that gives us a vector in Lambda^4 O
-- take one nonzero coordinate of this vector (preferable 
-- take a coordinate that is already 1)
-- localize equations to that coordinate patch
-- take partial derivatives of all the equations with respect
-- to all the local variables (i.e. compute the jacobian)
-- then "plug in the point"
-- see if the Jacobian is still of rank 1
-- funny thing is you probably won't need all the coordinates 
-- because some of the will be supplied by the plucker relations

wil0246 = indMatOnExtAlg (cols [0,2,4,6] changeOfBasis) 4
coorEqsOne = fst . head . filter ((==1) . snd) . zip allCoordinates $ head . transpose $ matrixAsList wil0246
-- it turns out coorEqsOne = "0246"
eqsAtProb = map simplify $ localizeDefEqs coorEqsOne
locVars' = localVars coorEqsOne
test v = if elem v locVars' then "T" else "FALSE"

isEverythingLocalized = map (mapVars test) eqsAtProb
-- THERE IS SOMETHING WRONG HERE!
-- CHECK WHAT'S UP WITH THIS!!!
--
-- ALSO, WE WERE USING A SHORTCUT WHEN COMPUTING JACOBIANS
-- CHECK THAT AS WELL! COMPARE IT TO ACTUALLY LOCALIZING
-- EQUATIONS AND THEN TAKING DERIVATIVES AND PLUGGING IN
-- ZEROES EVERYWHERE
