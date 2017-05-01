module LaurentPolynomials where

import Data.List
import Data.Maybe

type Variable = String
type Power = Int

-- a single term is a coefficient times the product of variables with powers
data LTerm a = Term a [(Variable, Power)]

-- a polynomial is sum of single terms as defined above
data LPoly a = Poly [LTerm a]

instance (Show a) => Show (LTerm a) where
  show (Term c ls) = concat ((show c):vars) where
    vars = map paste ls
    paste (var, 1) = ' ':var
    paste (var, pow) = ' ':var ++ "^{" ++ show pow ++ "}"

instance (Show a) => Show (LPoly a) where
  show (Poly []) = "0"
  show (Poly ls) = intercalate " + " $ map show ls

instance (Eq a, Num a) => Eq (LTerm a) where
  t1 == t2 = a == b && ls == ks where
    Term a ls = simplifyTerm t1
    Term b ks = simplifyTerm t2

instance (Eq a, Num a) => Eq (LPoly a) where
  p1 == p2 = ls1 == ls2 where
    Poly ls1 = simplify p1 
    Poly ls2 = simplify p2

instance (Eq a, Num a) => Num (LPoly a) where
  (Poly ls1) + (Poly ls2) = simplify $ Poly (ls1++ls2)
  p1 * p2 = simplify $ multLPolys p1 p2
  fromInteger n = Poly [Term (fromInteger n) []]
  negate (Poly ls) = Poly (map negTerm ls) where
    negTerm (Term c l) = Term (-c) l

instance Functor (LTerm) where
  fmap f (Term c ls) = Term (f c) ls

instance Functor (LPoly) where
  fmap f (Poly ls) = Poly (map (fmap f) ls)

instance (Num a, Ord a) => Ord (LTerm a) where
  compare (Term a1 ls1) (Term a2 ls2) = if cmp == EQ then compare a1 a2 else cmp where
    cmp = compare ls1 ls2

instance (Num a, Ord a) => Ord (LPoly a) where
  compare (Poly ls1) (Poly ls2) = compare ls1 ls2

simplifyTerm :: (Eq a, Num a) => LTerm a -> LTerm a
simplifyTerm (Term 0 _) = Term 0 []
simplifyTerm (Term a ls) = Term a newls where
  newls = filterOutZeroes . addExponents . groupBy varName $ sortOn fst ls
  varName (a,_) (b,_) = a == b
  addExponents = map g
  g ls = let (vars, pows) = unzip ls in (head vars, sum pows)
  filterOutZeroes = filter ((/=0) . snd)
  
simplify :: (Eq a, Num a) => LPoly a -> LPoly a
simplify (Poly ls) = Poly newls where
  newls = filterOutZeroes . addTerms . groupBy compVars . sortOn vars $ map simplifyTerm ls
  coeffs (Term a _) = a
  vars (Term _ ls) = ls
  compVars a b = (vars a) == (vars b)
  addTerms = map (\ls -> (Term (sum (map coeffs ls)) (vars $ head ls)))
  filterOutZeroes = filter ((/=0) . coeffs)

plugInValInTerm :: (Eq a, Num a) => String -> a -> LTerm a -> LTerm a
plugInValInTerm c v q   =  retVal where
  Term k ls             =  simplifyTerm q
  maybeIndex            =  findIndex ((==c) . fst) ls
  (init,(_,power):rest) =  splitAt (fromJust maybeIndex) ls
  retVal                =  if isNothing maybeIndex
                            then Term k ls
                            else Term (k*(v ^ power)) (init++rest)

plugInVal :: (Eq a, Num a) => String -> a -> LPoly a -> LPoly a
plugInVal c v (Poly ls) = simplify $ Poly (map (plugInValInTerm c v) ls)

multLTerms :: (Num a) => LTerm a -> LTerm a -> LTerm a
multLTerms (Term c1 ls1) (Term c2 ls2) = Term (c1*c2) (ls1++ls2)

multLTermLPoly :: (Num a) => LTerm a -> LPoly a -> LPoly a
multLTermLPoly t1 (Poly ls) = Poly $ map (multLTerms t1) ls

multLPolys :: (Num a) => LPoly a -> LPoly a -> LPoly a
multLPolys (Poly ls1) (Poly ls2) = Poly [multLTerms t1 t2 | t1 <- ls1, t2 <- ls2]

powerLPoly :: (Eq a, Num a) => Int -> LPoly a -> LPoly a
powerLPoly 0 _ = Poly [Term 1 []]
powerLPoly n p = (iterate (simplify . multLPolys p) p) !! (n-1)

plugInIdentityInTerm :: (Eq a, Num a) => String -> LPoly a -> LTerm a -> LPoly a
plugInIdentityInTerm c v q   =  retVal where
  Term k ls             =  simplifyTerm q
  maybeIndex            =  findIndex ((==c) . fst) ls
  (init,(_,power):rest) =  splitAt (fromJust maybeIndex) ls
  vpowered              =  powerLPoly power v
  retVal                =  if isNothing maybeIndex
                            then Poly [Term k ls]
                            else multLTermLPoly (Term k (init++rest)) vpowered

polyAsList :: LPoly a -> [LTerm a]
polyAsList (Poly ls) = ls

sumPolys :: [LPoly a] -> LPoly a
sumPolys [] = Poly []
sumPolys (Poly ls:rest) = Poly $ ls ++ (polyAsList $ sumPolys rest)

plugInIdentity :: (Eq a, Num a) => String -> LPoly a -> LPoly a -> LPoly a
plugInIdentity c v (Poly ls) = sumPolys $ map (plugInIdentityInTerm c v) ls

mapVarsInTerms :: (String -> String) -> LTerm a -> LTerm a
mapVarsInTerms f (Term c ls) = Term c (map (\(v,p) -> (f v,p) ) ls)

mapVars :: (String -> String) -> LPoly a -> LPoly a
mapVars f (Poly ls) = Poly $ map (mapVarsInTerms f) ls

