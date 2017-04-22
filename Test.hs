
import ComplexNumbers (ComplexNumber, i)
import ComplexOctonions
import CayleyPlanes
import DifferentialFormCalculus

nulls = [1+ii*oi, oi+ii*oj, oj + ii*ok]

isNormEqTrueForQuadCrPr = map (normSquared . quadCrPrOnList) gen4VsIn8d == replicate (length gen4VsIn8d) 1

capPhiSquared = wdgPrOnComps capPhiComps capPhiComps

sumOfCapPhiAndImQd = map (\q -> ((capPhiOnList q)**2) + normSquared (imaginary $ preQuadCrProdOnList q)) gen4VsIn8d

isAllOnes :: [ComplexNumber] -> Bool
isAllOnes [] = True
isAllOnes (1:xs) = isAllOnes xs
isAllOnes _ = False

