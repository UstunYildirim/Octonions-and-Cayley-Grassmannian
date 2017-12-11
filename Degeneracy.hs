module Degeneracy where

import Data.List
import ComplexNumbers (ComplexNumber, i)
import ComplexOctonions
import LaurentPolynomials
import DifferentialFormCalculus
import MatrixCalculus
import CayleyPlanes
import MaximalTorus
import Plucker
import Jacobians


itimes = scalarMult i
newBasis = [1 + itimes oi, 1 - itimes oi, 
            oj + itimes ok, oj - itimes ok,
            ol + itimes oli, ol - itimes oli,
            olj + itimes olk, olj - itimes olk]

choose n l = uniqStrIncNtuple n $ listAsNewBasisVecs l

listAsNewBasisVecs l = sublist (map (\c-> read [c] :: Int) l) newBasis

matRepofB l = let l' = listAsNewBasisVecs l in
  [[dotProd u v | u <- l'] | v <- l']

bSing = map matRepofB singularFxdPnts
bReg = map matRepofB regularFxdPnts
tripSing = map (map triCrProdOnList . choose 3) singularFxdPnts
tripReg = map (map triCrProdOnList . choose 3) regularFxdPnts

countZeroRows [] = 0
countZeroRows (r:rs) = if r == [0,0,0,0] then 1 + countZeroRows rs else countZeroRows rs
