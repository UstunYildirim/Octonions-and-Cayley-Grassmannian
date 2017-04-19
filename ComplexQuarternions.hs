module ComplexQuarternions where

import ComplexNumbers hiding (coordinateWiseOperation, sumCoordinates, multiply, normSquared, conj, add, norm, scalarMult, real)
import qualified Data.List as List

q1 :: Quarternion
q1 = Quarternion(1,0,0,0)
qi :: Quarternion
qi = Quarternion(0,1,0,0)
qj :: Quarternion
qj = Quarternion(0,0,1,0)
qk :: Quarternion
qk = Quarternion(0,0,0,1)

data Quarternion = Quarternion (ComplexNumber,ComplexNumber,ComplexNumber,ComplexNumber)
                                  deriving (Eq)
                                  
quarternionAsList :: Quarternion -> [ComplexNumber]
quarternionAsList (Quarternion (a,b,c,d)) = [a,b,c,d]

coordinateWiseOperation :: (ComplexNumber -> ComplexNumber -> ComplexNumber) -> Quarternion -> Quarternion -> Quarternion
coordinateWiseOperation f (Quarternion(a,b,c,d)) (Quarternion(x,y,z,w))
  = Quarternion(f a x, f b y, f c z, f d w)
add :: Quarternion -> Quarternion -> Quarternion
add = coordinateWiseOperation (+)

subtract :: Quarternion -> Quarternion -> Quarternion
subtract = coordinateWiseOperation (-)

-- this is coordinatewise!
multiply :: Quarternion -> Quarternion -> Quarternion
multiply = coordinateWiseOperation (*)

sumCoordinates :: Quarternion -> ComplexNumber
sumCoordinates (Quarternion (a,b,c,d)) = a+b+c+d

real :: Quarternion -> ComplexNumber
real (Quarternion (a,_,_,_)) = a

imaginary :: Quarternion -> Quarternion
imaginary (Quarternion (_,b,c,d)) = Quarternion(0,b,c,d)

isReal :: Quarternion -> Bool
isReal (Quarternion (_,0,0,0)) = True
isReal _ = False

conj :: Quarternion -> Quarternion
conj (Quarternion (a,b,c,d)) = (Quarternion (a,-b,-c,-d))

dotProd :: Quarternion -> Quarternion -> ComplexNumber
dotProd q q' = real $ q * (conj q')

normSquared :: Quarternion -> ComplexNumber
normSquared (Quarternion (a,b,c,d)) = a*a + b*b + c*c + d*d -- dotProd q q

norm :: Quarternion -> ComplexNumber
norm = sqrt . normSquared

scalarMult :: ComplexNumber -> Quarternion -> Quarternion
scalarMult s (Quarternion (a,b,c,d)) = Quarternion(s*a,s*b,s*c,s*d)

qmult :: Quarternion -> Quarternion -> Quarternion
qmult (Quarternion (a,b,c,d)) (Quarternion (x,y,z,w)) =
  let
    a' = a*x-b*y-c*z-d*w
    b' = a*y+b*x+c*w-d*z
    c' = a*z-b*w+c*x+d*y
    d' = a*w+b*z-c*y+d*x
  in
    Quarternion (a',b',c',d')

instance Num Quarternion where
 (+) = add
 (*) = qmult
 (-) = ComplexQuarternions.subtract
 negate (Quarternion q1) = scalarMult (-1) (Quarternion q1)
 fromInteger c = (Quarternion (fromInteger c,0,0,0))
 abs q = Quarternion (norm q,0,0,0)
 signum (Quarternion (a,_,_,_)) = Quarternion (signum a,0,0,0)

instance Fractional Quarternion where
  fromRational r = Quarternion (fromRational r,0,0,0)
  recip q@(Quarternion (a,b,c,d)) = Quarternion(a',b',c',d') where
    n = normSquared q
    a' = a/n
    b' = -b/n 
    c' = -c/n 
    d' = -d/n 

instance Show Quarternion where
  show q = List.intercalate " + " $ map (\(x,y)->x++y) $ filter (\(x,_)->x/="(0)") $ zip (map (\x->'(':(show x)++")") (quarternionAsList q)) ["","i","j","k"] 

instance Ord Quarternion where
  q < q' = (quarternionAsList q) < (quarternionAsList q')
  q <= q' = (quarternionAsList q) <= (quarternionAsList q')

instance Floating Quarternion where
  q ** n 
    | n == 0 = Quarternion(1,0,0,0)
    | n > 0 = q*(q**(n-1))
    | otherwise = recip (q ** (-n))

