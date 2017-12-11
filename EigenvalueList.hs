module EigenvalueList where

import Data.List
import ComplexNumbers (ComplexNumber, i)
import ComplexOctonions
import LaurentPolynomials
import DifferentialFormCalculus
import MatrixCalculus
import CayleyPlanes
import MaximalTorus
import Plucker

la :: LPoly (ComplexNumber)
la  = Poly [Term 1 [("\\lambda",1)]]
lai = Poly [Term 1 [("\\lambda",-1)]]
mu  = Poly [Term 1 [("\\mu",1)]]
mui = Poly [Term 1 [("\\mu",-1)]]
ga  = Poly [Term 1 [("\\gamma",1)]]
gai = Poly [Term 1 [("\\gamma",-1)]]

eigvals = [la*mu, lai*mui,
           la*mui, lai*mu,
           la*ga, lai*gai,
           la*gai, lai*ga]

indEigVals = map product $ uniqStrIncNtuple 4 eigvals

eigValLookUp var = lookup var pairWith4Vecs'

pairWith4Vecs' = zip allCoordinates indEigVals
pairWith4Vecs = zip indEigVals allCoordinates

vecWithTilde = ("\\wil e_{" ++) . (++ "}")

collectSameEigVals :: [(LPoly ComplexNumber, String)] -> [(LPoly ComplexNumber, [String])]
collectSameEigVals = map createEigValVectorListTuple . groupBy fsts . sortOn fst
  where
    fsts a b = (fst a) == (fst b)
    createEigValVectorListTuple ls = 
      let (eigval:_, vs) = unzip ls in (eigval, vs)

latexTableOfEigValsAndVecs pairs = intercalate "\\\\ \\hline\n" $ map genLine (collectSameEigVals pairs)
  where
    genLine (val, vecs) = let vecs' = intercalate ", " $ map vecWithTilde vecs in 
      " $" ++ show val ++ "$ & $" ++ vecs' ++ "$ "

-- OLD
-- eigVecsInXminWithVals = map (\(n,vec) -> (indEigVals !! n, vec)) . filter test $ zip [0..] allCoordinates
--   where
--     test (n,_) = col n newCoordsImQdCrPr == zeroMatrix 8 1

itimes = scalarMult i
newBasis = [1 + itimes oi, 1 - itimes oi, 
            oj + itimes ok, oj - itimes ok,
            ol + itimes oli, ol - itimes oli,
            olj + itimes olk, olj - itimes olk]
eigVecsInXminWithVals = map (\(n,vec) -> (indEigVals !! n, allCoordinates !! n)) . filter ((==0) .normSquared . imQuadCrPrOnList . snd) . zip [0..] $ ngenFromBasis 4 newBasis

eigVecsInXmin = map snd eigVecsInXminWithVals
