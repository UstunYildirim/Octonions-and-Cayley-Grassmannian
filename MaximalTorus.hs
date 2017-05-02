module MaximalTorus where

import ComplexNumbers (ComplexNumber, i)
import ComplexOctonions
import LaurentPolynomials
import DifferentialFormCalculus
import MatrixCalculus
import CayleyPlanes
import Plucker

simplifyMatrix :: (Eq a, Num a) => Matrix (LPoly a) -> Matrix (LPoly a)
simplifyMatrix = Matrix . map (map simplify) . matrixAsList

capPhiAsRowMatrix = Matrix [map (value . flip lookup capPhiComps) basis] -- in the basis e^{ijkl} dual to e_{ijkl}
  where
    basis = allCoordinates
    value Nothing = Poly []
    value (Just a) = Poly [Term a []]


doesInvMatPreserveCapPhi m = (capPhiAsRowMatrix * indMatOnExtAlg m 4) == capPhiAsRowMatrix

p var = Poly [Term (1/2) [(var,1)], Term (1/2) [(var,-1)]]
m var = Poly [Term (1/2) [(var,1)], Term (-1/2) [(var,-1)]]
half = Poly [Term (0.5) []]
quarter = Poly [Term (0.25) []]
halfQuarter = Poly [Term (0.125) []]
ipoly = Poly [Term (i) []]
halfiPoly = Poly [Term (0.5*i) []]

lMat :: String -> Matrix (LPoly ComplexNumber)
lMat var = Matrix [[p var, (-ipoly)*m var],
                   [ipoly*m var, p var]]

lMatInv :: String -> Matrix (LPoly ComplexNumber)
lMatInv var = Matrix [[p var, ipoly*m var],
                   [(-ipoly)*m var, p var]]

ll = lMat "lambda"
lli = lMatInv "lambda"
lm = lMat "mu"
lmi = lMatInv "mu"
lg = lMat "gamma"
lgi = lMatInv "gamma"

aMat = matDirSum m m where
  m = matDirSum ll ll
bMat = matDirSum m (identity 4) where
  m = matDirSum lm lmi
cMat = matDirSum (identity 4) m where
  m = matDirSum lg lgi

doesInvOfApreserveCapPhi = doesInvMatPreserveCapPhi aMat
doesInvOfBpreserveCapPhi = doesInvMatPreserveCapPhi bMat
doesInvOfCpreserveCapPhi = doesInvMatPreserveCapPhi cMat

chbas = Matrix [[1,     1],
                [ipoly, -ipoly]]
chbasInv = Matrix [[half, -halfiPoly],
                   [half, halfiPoly]]
chbasOp = Matrix [[1,      1],
                  [-ipoly, ipoly]]

changeOfBasis = matDirSum (matDirSum chbas chbas) (matDirSum chbas chbas)

-- the following is an 7x70 matrix representing the imaginary part of the 
-- quadruple cross product. 
imQdCrPrAsMat = Matrix [map (value . flip lookup (pickCoeffsIQCP v)) basis | v <- octGens]
  where
    basis = map (("e^{" ++) . (++ "}")) allCoordinates
    value Nothing = Poly []
    value (Just a) = Poly [Term a []]


newCoordsImQdCrPr = imQdCrPrAsMat * (indMatOnExtAlg changeOfBasis 4)

transformedDefEqs = map (filter (\(_,a)-> a/=0)) . map (zip allCoordinates) $ matrixAsList newCoordsImQdCrPr

normTransDefEqs = zipWith zipper [-halfQuarter*ipoly, quarter, -quarter*ipoly, quarter, -quarter*ipoly, -quarter, quarter*ipoly] transformedDefEqs
  where
    zipper scalar = map (\(var, coeff) -> (var, scalar * coeff))


latexPrintRescaledEqs = map latexPrintEq normTransDefEqs
  where
    latexPrintEq = concatMap processVarAndCoeff
    processVarAndCoeff (var, 1) = " + \\wil p_{" ++ var ++ "}"
    processVarAndCoeff (var, -1) = " - \\wil p_{" ++ var ++ "}"
    processVarAndCoeff _ = error "coefficients other than +1, -1 are not handled!"
