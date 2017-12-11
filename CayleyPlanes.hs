module CayleyPlanes where

import ComplexOctonions
import ComplexNumbers (ComplexNumber, i)
import Data.List
import DifferentialFormCalculus
import LaurentPolynomials

gen3VsIn7d :: [[Octonion]]
gen3VsIn8d :: [[Octonion]]
gen4VsIn7d :: [[Octonion]]
gen4VsIn8d :: [[Octonion]]
gen3VsIn7d = ngenFromBasis 3 imagOctGens
gen3VsIn8d = ngenFromBasis 3 octGens
gen4VsIn7d = ngenFromBasis 4 imagOctGens
gen4VsIn8d = ngenFromBasis 4 octGens

vName :: Octonion -> Char
vName = unjust . flip lookup (zip octGens "01234567")
  where
    unjust (Just v) = v
    unjust Nothing  = error "non generator vector being named!"

vNames :: [Octonion] -> String
vNames = map vName

nvecNames :: [Octonion] -> String
ncovNames :: [Octonion] -> String
nvecNames os = "e_{" ++ vNames os ++ "}"
ncovNames os = "e^{" ++ vNames os ++ "}"

removeVecCov = init . drop 3

formComps7 n f = map (\(a,b) -> (vNames a, b)) $ formComps n imagOctGens f
formComps8 n f = map (\(a,b) -> (vNames a, b)) $ formComps n octGens f
formComps7' n f = map (\(a,b) -> (ncovNames a, b)) $ formComps n imagOctGens f
formComps8' n f = map (\(a,b) -> (ncovNames a, b)) $ formComps n octGens f

preCrProdOnList :: [Octonion] -> Octonion
preCrProdOnList = alternate f where
  f [u, v] = (conj v) * u

crProd :: Octonion -> Octonion -> Octonion
crProd u v = imaginary ((conj v) * u)

crProdOnList :: [Octonion] -> Octonion
crProdOnList [u,v] = crProd u v

phi :: Octonion -> Octonion -> Octonion -> ComplexNumber
phi u v w = dotProd u (crProd v w)

phiOnList :: [Octonion] -> ComplexNumber
phiOnList [u,v,w] = phi u v w
phiOnList _ = error "phi requires exactly 3 vectors"

phiComps = formComps7 3 phiOnList

contractPhi x = formComps7 2 (intPr x phiOnList)

metricAndVol x = wdgPrOnComps (contractPhi x) $ wdgPrOnComps (contractPhi x) phiComps

associator :: Octonion -> Octonion -> Octonion -> Octonion
associator u v w = scalarMult 0.5 $ u * (v * w) - (u * v) * w

associatorOnList :: [Octonion] -> Octonion
associatorOnList [u,v,w] = associator u v w

-- we want to express triCrProd as alternation
-- of u(\bar v w).

preTrCrPrOnList = alternate f where
  f [u, v, w] = (u * (conj v)) * w

-- this definition of triple cross product agrees with the one below.

triCrProd :: Octonion -> Octonion -> Octonion -> Octonion
triCrProd u v w = (scalarMult 0.5) ((u * cnjv) * w - (w * cnjv) * u) where cnjv = conj v

triCrProdOnList :: [Octonion] -> Octonion
triCrProdOnList [u,v,w] = triCrProd u v w
triCrProdOnList _ = error "triCrProd requires exactly 3 vectors"

capPhi :: Octonion -> Octonion -> Octonion -> Octonion -> ComplexNumber
capPhi x u v w = dotProd x (triCrProd u v w)

capPhiOnList :: [Octonion] -> ComplexNumber
capPhiOnList [x,u,v,w] = capPhi x u v w
capPhiOnList _ = error "capPhi is a 4 form! It requires 4 vectors to be evaluated."

capPhiComps = formComps8 4 capPhiOnList

latexPrintCapPhi = latexPrintFormNum $ formComps8' 4 capPhiOnList

preQuadCrProdOnList = alternate f
  where
    --f [x,u,v,w] = negate $ ((x * conj u) * v) * conj w
    f [x,u,v,w] = - triCrProdOnList [x,u,v] * conj w

quadCrPr x u v w = scalarMult (-0.25) (
  (triCrProd x u v) * conj w
  -(triCrProd w x u) * conj v
  +(triCrProd v w x) * conj u
  -(triCrProd u v w) * conj x)

quadCrPrOnList [x,u,v,w] = quadCrPr x u v w

imQuadCrPrOnList = imaginary . quadCrPrOnList

imQdCrPrComps = formComps8' 4 imQuadCrPrOnList

pickCoeffsIQCP o = map (\(x,y)->(x,real (y/o))) $ filter (isMultiple o . snd) imQdCrPrComps

latexPrintImQdCrPrComp o = (++"\\right)"++ (nvecNames [o])) . ("\\left("++) . latexPrintFormNum $ pickCoeffsIQCP o

definingEquation o = Poly . map createTerm $ pickCoeffsIQCP o where
    createTerm (x,y) = Term (y) [(x',1)]
      where x' = init $ drop 3 x

definingEquations = map definingEquation imagOctGens

