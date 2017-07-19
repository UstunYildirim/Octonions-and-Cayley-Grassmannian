module Check where

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

allOct2 = [[u,v] | u<-octGens, v<-octGens]
allOct3 = [[u,v,w] | u<-octGens, v<-octGens, w<-octGens]
allOct4 = [[x,u,v,w] | x<-octGens, u<-octGens, v<-octGens, w<-octGens]
wedgeOct2 = uniqStrIncNtuple 2 octGens
wedgeOct3 = uniqStrIncNtuple 3 octGens
wedgeOct4 = uniqStrIncNtuple 4 octGens

doConventionsAgreeWithPaper = all (==True)
  [ isCompositionAlgebra,
    map normSquared octGens == replicate 8 1,
    map (^2) octGens == 1:replicate 7 (-1),
    map (*1) octGens == octGens,
    oi*oj == ok,
    oj*ok == oi,
    ok*oi == oj,
    ol*oi == oli,
    ol*oj == olj,
    ol*ok == olk,
    olj*oi == olk,
    olk*oj == oli,
    oli*ok == olj,
    map crProdOnList allOct2 == map preCrProdOnList allOct2,
    map triCrProdOnList allOct3 == map preTrCrPrOnList allOct3,
    map quadCrPrOnList allOct4 == map preQuadCrProdOnList allOct4
  ] -- this is true

doConventionsAgreeWithPaper2 = all (==True)
  [ map (normSquared . quadCrPrOnList) wedgeOct4 == replicate 70 1,
    map phiOnList wedgeOct3 == map (real . triCrProdOnList) wedgeOct3,
    map capPhiOnList wedgeOct4 == map (real . quadCrPrOnList) wedgeOct4,
    -- phi'nin power'i dogru volume formu veriyor mu?
    -- capPhi paper'da yazan sekilde mi gozukuyor?
    doesInvOfApreserveCapPhi,
    doesInvOfBpreserveCapPhi,
    doesInvOfCpreserveCapPhi
  ]
