module ComplexOctonions where

import ComplexNumbers hiding (coordinateWiseOperation, sumCoordinates, multiply, normSquared, conj, add, norm, scalarMult, real, imaginary, isReal)
import ComplexQuarternions hiding (coordinateWiseOperation, sumCoordinates, multiply, normSquared, conj, add, norm, scalarMult, real,dotProd, imaginary, isReal)
import qualified ComplexQuarternions as CQ
import qualified Data.List as List

ii  :: Octonion 
oi  :: Octonion 
oj  :: Octonion 
ok  :: Octonion 
ol  :: Octonion 
oli :: Octonion 
olj :: Octonion 
olk :: Octonion 
ii  = scalarMult i 1
oi  = Octonion(qi,0)
oj  = Octonion(qj,0)
ok  = Octonion(qk,0)
ol  = Octonion(0,1)
oli = ol * oi
olj = ol * oj
olk = ol * ok

octGens     :: [Octonion] 
ioctGens    :: [Octonion] 
imagOctGens :: [Octonion] 
octGens     = [1,oi,oj,ok,ol,oli,olj,olk]
ioctGens    = map (scalarMult i) octGens
imagOctGens = tail octGens

data Octonion = Octonion {  octonion :: (Quarternion,Quarternion)
                                  } deriving (Eq)
octonionAsList :: Octonion -> [Quarternion]
octonionAsList (Octonion (a,b)) = [a,b]

octonionAsCxList :: Octonion -> [ComplexNumber]
octonionAsCxList (Octonion (Quarternion (a,b,c,d),Quarternion(x,y,z,w))) = [a,b,c,d,x,-y,-z,-w]

coordinateWiseOperation :: (Quarternion -> Quarternion -> Quarternion) -> Octonion -> Octonion -> Octonion
coordinateWiseOperation f (Octonion(a,b)) (Octonion(x,y)) = Octonion(f a x, f b y)

add :: Octonion -> Octonion -> Octonion
subtract :: Octonion -> Octonion -> Octonion
multiply :: Octonion -> Octonion -> Octonion

add = coordinateWiseOperation (+)
subtract = coordinateWiseOperation (-)
multiply = coordinateWiseOperation (CQ.multiply)

scalarMult :: ComplexNumber -> Octonion -> Octonion
scalarMult c (Octonion(a,b)) = Octonion(qsm c a, qsm c b) where
  qsm = CQ.scalarMult

sumCoordinates :: Octonion -> ComplexNumber
sumCoordinates (Octonion (a,b)) = (qsc a) + (qsc b)
  where
    qsc = CQ.sumCoordinates

conj :: Octonion -> Octonion
conj (Octonion (a,b)) = Octonion (CQ.conj a, -b)

normSquared :: Octonion -> ComplexNumber
omult :: Octonion -> Octonion -> Octonion

normSquared (Octonion (Quarternion (a,b,c,d),Quarternion(x,y,z,w))) = a*a + b*b + c*c + d*d + x*x + y*y + z*z + w*w 
omult (Octonion (x,y)) (Octonion (u,v)) = Octonion ((x`qmult`u) `CQ.subtract` ((CQ.conj v)`qmult`y),
                                                    (v`qmult`x) `CQ.add` (y`qmult`(CQ.conj u)))

norm :: Octonion -> ComplexNumber
norm = sqrt . normSquared

real :: Octonion -> ComplexNumber
imaginary :: Octonion -> Octonion

real (Octonion (a,_)) = CQ.real a
imaginary (Octonion (a,b)) = Octonion(CQ.imaginary a, b)

dotProd :: Octonion -> Octonion -> ComplexNumber
dotProd u v = real (u * conj v)

dotProdOnList [u,v] = dotProd u v

isMultiple o1 o2 = isReal (o2/o1)

instance Num Octonion where
  (+)                     = add
  (*)                     = omult
  (-)                     = ComplexOctonions.subtract
  negate (Octonion q1)    = scalarMult (-1) (Octonion q1)
  fromInteger c           = (Octonion (fromInteger c,0))
  abs q                   = Octonion ((CQ.scalarMult (normSquared q) 1),0)
  signum (Octonion (a,_)) = Octonion (signum a,0)

fromComplexNumber :: ComplexNumber -> Octonion
fromComplexNumber cn = Octonion (Quarternion (cn,0,0,0), 0)

instance Fractional Octonion where
  fromRational r = Octonion (fromRational r,0)
  recip q = scalarMult (1/(normSquared q)) (conj q)

instance Show Octonion where
  show 0 = "0"
  show q = List.intercalate " + " .
    map (\(x,y)->y++x) .
      filter (\(x,_)->x/="[]") $
        zip ([a', b'])  ["","l"] 
         where
          [a,b] = octonionAsList q
          a' = putBrackets $ show a
          b' = putBrackets $ show (CQ.conj b)
          putBrackets = ('[':) . (++"]")

instance Ord Octonion where
  q < q' = (octonionAsList q) < (octonionAsList q')
  q <= q' = (octonionAsList q) <= (octonionAsList q')

instance Floating Octonion where
  q ** n 
    | n == 0 = Octonion(1,0)
    | n > 0 = q*(q**(n-1))
    | otherwise = recip (q ** (-n))

isCompositionAlgebra = and [conj (u*v) == (conj v)*(conj u)| u<-octGens, v<-octGens]

isReal :: Octonion -> Bool
isReal (Octonion (a,0)) = CQ.isReal a
isReal _ = False

makeUnit :: Octonion -> Octonion -- in the sense that norm = 1 note: there are null vectors!
makeUnit q = scalarMult (1/(norm q)) q
