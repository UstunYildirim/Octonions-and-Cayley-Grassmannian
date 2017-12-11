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
