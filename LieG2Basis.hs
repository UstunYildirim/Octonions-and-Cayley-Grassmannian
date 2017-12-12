module LieG2Basis where

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


twiceIdentity = Matrix [[if (i==j) then 2 else 0 | j <- [0..7] ] | i <- [0..7]] :: Matrix ComplexNumber

elemMatrix i0 j0 = Matrix [[if ((i0 == i) && (j0 == j)) then 1 else 0 | j <- [0..7]] | i <- [0..7]]
skewSymMat i0 j0 = elemMatrix i0 j0 - elemMatrix j0 i0
l = skewSymMat

chBasBlk = Matrix [[1, 1], [i, -i]]
chBasBlkInv = Matrix [[1/2, -i/2], [1/2, i/2]]
chBas = d a a where a = d b b
                    d = matDirSum
                    b = chBasBlk
chBasInv = d a a where a = d b b
                       d = matDirSum
                       b = chBasBlkInv

g2Basis = [l 1 6 + l 2 5, l 2 5 + l 3 4,
           l 1 7 + l 2 4, l 1 7 + l 3 5,
           l 1 4 + l 3 6, l 2 7 + l 3 6,
           l 1 5 - l 2 6, l 1 5 - l 3 7,
           l 1 2 + l 4 7, l 1 2 + l 5 6,
           l 1 3 + l 5 7, l 4 6 + l 5 7,
           l 2 3 + l 4 5, l 2 3 + l 6 7 ]

g2TildeBasis = map (\x -> chBasInv * x * chBas) g2Basis

cxToSymWvar v c = Poly [Term c [([v],1)]]

symg2Basis = map (\(c,m) -> simplifyMatrix $ fmap (cxToSymWvar c) m) $ zip ['a'..'z'] g2Basis
singleMatStandardBasis = foldr1 (+) symg2Basis

x i = g2Basis !! (i-1)
xt i = g2TildeBasis !! (i-1)

divi = matScalarMult (1/i)
-- generators of the stabilizer subalgebra
lieParabolic = [divi (x 1) + x 4, divi (x 1) + x 3 - divi (x 2),
                divi (x 5) + x 7, divi (x 5) + x 8 - divi (x 6),
                divi (x 9) + x 11, divi (x 10) + x 11,
                x 12, divi (x 13), divi (x 14)]
lieParabolicTilde = [divi (xt 1) + xt 4, divi (xt 1) + xt 3 - divi (xt 2),
                divi (xt 5) + xt 7, divi (xt 5) + xt 8 - divi (xt 6),
                divi (xt 9) + xt 11, divi (xt 10) + xt 11,
                xt 12, divi (xt 13), divi (xt 14)]
symLieParabolic = map (\(c,m) -> simplifyMatrix $ fmap (cxToSymWvar c) m) $ zip ['A'..'Z'] lieParabolic
symLieParabolicTilde = map (\(c,m) -> simplifyMatrix $ fmap (cxToSymWvar c) m) $ zip ['a'..'z'] lieParabolicTilde
singleMatParabolicStandardBasis = foldr1 (+) symLieParabolic
singleMatParabolicTildeBasis = foldr1 (+) symLieParabolicTilde

