module ComplexNumbers where

import qualified Data.List as List

i :: ComplexNumber
i = ComplexNumber(0,1)

data ComplexNumber = ComplexNumber (Double,Double) deriving (Eq)

complexNumberAsList :: ComplexNumber -> [Double]
complexNumberAsList (ComplexNumber (a,b)) = [a,b]

coordinateWiseOperation :: (Double -> Double -> Double) -> ComplexNumber -> ComplexNumber -> ComplexNumber
coordinateWiseOperation f (ComplexNumber (x1,y1)) (ComplexNumber (x2,y2)) = ComplexNumber (f x1 x2, f y1 y2) 

add :: ComplexNumber -> ComplexNumber -> ComplexNumber
add = coordinateWiseOperation (+)

subtract :: ComplexNumber -> ComplexNumber -> ComplexNumber
subtract = coordinateWiseOperation (-)

-- this is coordinatewise!
multiply :: ComplexNumber -> ComplexNumber -> ComplexNumber
multiply = coordinateWiseOperation (*)

sumCoordinates :: ComplexNumber -> Double
sumCoordinates (ComplexNumber (a,b)) = a+b

normSquared :: ComplexNumber -> Double
normSquared q = real $ q * (conj q)

norm :: ComplexNumber -> Double
norm = sqrt . normSquared

scalarMult :: Double -> ComplexNumber -> ComplexNumber
scalarMult c (ComplexNumber (a,b)) = ComplexNumber (c*a,c*b)

conj :: ComplexNumber -> ComplexNumber
conj (ComplexNumber (x,y)) = ComplexNumber (x,-y)

cmult :: ComplexNumber -> ComplexNumber -> ComplexNumber
cmult (ComplexNumber (a,b)) (ComplexNumber (x,y)) = ComplexNumber (a*x-b*y,a*y+b*x)

real :: ComplexNumber -> Double
real (ComplexNumber (a,_)) = a

imaginary :: ComplexNumber -> Double
imaginary (ComplexNumber (_,b)) = b

isReal :: ComplexNumber -> Bool
isReal (ComplexNumber (_,0)) = True
isReal _ = False

isImaginary :: ComplexNumber -> Bool
isImaginary (ComplexNumber (0,_)) = True
isImaginary _ = False

instance Num ComplexNumber where
 (+) = add
 (*) = cmult
 (-) = ComplexNumbers.subtract
 negate (ComplexNumber q1) = scalarMult (-1) (ComplexNumber q1)
 fromInteger c = (ComplexNumber (fromInteger c,0))
 abs q = ComplexNumber (norm q,0)
 signum (ComplexNumber (a,_)) = ComplexNumber (signum a,0)

instance Fractional ComplexNumber where
  fromRational r = ComplexNumber (fromRational r,0)
  recip q@(ComplexNumber (a,b)) = ComplexNumber(a',b') where
    n = normSquared q
    a' = a/n
    b' = -b/n 

instance Show ComplexNumber where
  show (ComplexNumber (0,0)) = "0"
  show (ComplexNumber (a,0)) = show a
  show (ComplexNumber (0,b)) = show b ++ "i"
  show (ComplexNumber (a,b)) = show a ++ " + " ++ show b ++ "i" 

instance Ord ComplexNumber where
  q < q' = (complexNumberAsList q) < (complexNumberAsList q')
  q <= q' = (complexNumberAsList q) <= (complexNumberAsList q')

instance Floating ComplexNumber where
  q ** n 
    | n == 0 = ComplexNumber (1,0)
    | n > 0 = q*(q**(n-1))
    | otherwise = recip (q ** (-n))

  sqrt 0 = 0
  sqrt q = let
      normq = (norm q)
      sqrtNormq = sqrt normq
      unitq = scalarMult (1 / normq) q
      findAngle q'@(ComplexNumber (a,b)) = if b < 0 then pi + (findAngle (-q')) else acos a
      halfAngle = (findAngle unitq)/2
    in
      scalarMult sqrtNormq (ComplexNumber (cos halfAngle, sin halfAngle))


