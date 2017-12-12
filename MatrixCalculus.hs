module MatrixCalculus where

import Data.List
import DifferentialFormCalculus

data Matrix a = Matrix [[a]] deriving (Eq)

instance (Show a) => Show (Matrix a) where
  show = intercalate "\n" . map ( ('[':) . (++ "]") . intercalate ",\t" . map show). matrixAsList

instance Functor (Matrix) where
  fmap f = Matrix . map (map f) . matrixAsList

instance (Eq a, Num a) => Num (Matrix a) where
  (Matrix m1) + (Matrix m2) = Matrix $ zipWith (zipWith (+)) m1 m2
  q1@(Matrix m1) * q2@(Matrix m2) = Matrix [[entryOfProdMat i j | j <- [0..nc-1] ] | i <- [0..nr-1] ]
    where
      entryOfProdMat i j = sum $ zipWith (*) ithRow jthCol where
        ithRow = m1 !! i
        jthCol = (transpose m2) !! j
      nr = numRow q1
      nc = numCol q2
  fromInteger n = Matrix [[fromInteger n]]
  negate = fmap negate

matScalarMult l m = fmap (l *) m

matrixAsList (Matrix ls) = ls
listAsMatrix ls = Matrix ls

transposeMat :: Matrix a -> Matrix a
transposeMat (Matrix m) = Matrix (transpose m)

numRow = length . matrixAsList
numCol = length . head . matrixAsList

entry :: Int -> Int -> Matrix a -> a
entry i j = (!!j) . (!!i) . matrixAsList

zeroMatrix :: (Num a) => Int -> Int -> Matrix a
zeroMatrix n m = Matrix . replicate n $ replicate m 0

identity :: (Num a) => Int -> Matrix a
identity n = Matrix [[if i == j then 1 else 0 | j <- [1..n] ] | i <- [1..n] ]

sublist :: [Int] -> [a] -> [a]
sublist indexSet = map snd . filter (\(i,_) -> elem i indexSet) . zip [0..]

submat :: [Int] -> [Int] -> Matrix a -> Matrix a
submat i j = Matrix . map (sublist j) . sublist i . matrixAsList

rows i = submat i [0..]
cols j = submat [0..] j

row i = rows [i]
col j = cols [j]

trace :: (Num a) => Matrix a -> a
trace (Matrix m) = sum $ zipWith (\r i-> r!!i) m [0..]

sub :: [a] -> [[a]]
sub [x] = [[]]
sub (x:xs) = xs:(map (x:) (sub xs))

submats :: [[a]] -> [[[[a]]]]
submats [[x]] = [[[[]]]]
submats x = [transpose (map sub subm) | subm<-(sub x)]

det :: (Eq a, Num a) => Matrix a -> a
det (Matrix [[x]]) = x
det (Matrix x) = sum (zipWith clevMult (head x) (zipWith (*) detmatrix sign))
  where
    -- clevMult (obviously) avoids computation of second term
    -- if the first term is already zero. this is particularly
    -- good when dealing with large (or large number of) sparse
    -- matrices as it avoids computation of determinant of
    -- (n-1)x(n-1) matrices
    clevMult 0 _ = 0
    clevMult a b = a * b
    detmatrix = [det (Matrix subm) | subm<-(head (submats x))]
    sign = cycle [1, (-1)]


-- the following takes matrices a and b --\\     a 0
-- and creates the block matrix         --//     0 b
matDirSum :: (Num a) => Matrix a -> Matrix a -> Matrix a
matDirSum qa@(Matrix a) qb@(Matrix b) = Matrix $ zipWith (++) leftSide rightSide where
  leftSide = a ++ bottomLeft
  rightSide = topRight ++ b
  topRight = matrixAsList $ zeroMatrix nra ncb
  bottomLeft = matrixAsList $ zeroMatrix nrb nca
  nra = numRow qa
  nca = numCol qa
  nrb = numRow qb
  ncb = numCol qb

minor :: (Eq a, Num a) => [Int] -> [Int] -> Matrix a -> a
minor i j = det . submat i j

indMatOnExtAlg :: (Eq a, Num a) => Matrix a -> Int -> Matrix a
indMatOnExtAlg _ 0 = 1
indMatOnExtAlg m n = Matrix [[minor i j m | j <- lsCol ] | i <- lsRow ] where
  lsRow = uniqStrIncNtuple n [0.. numRow m-1]
  lsCol = uniqStrIncNtuple n [0.. numCol m-1]

latexPrintMatrix :: (Show a) => Matrix a -> String
latexPrintMatrix = ("\\begin{pmatrix}\n" ++) .(++ "\n\\end{pmatrix}\n") . intercalate " \\\\ \n" . map (intercalate " & " . map show) . matrixAsList
