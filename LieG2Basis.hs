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
