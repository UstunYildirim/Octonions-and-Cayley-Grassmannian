-- this module is currently useless!
module SingularTest where

import System.IO
import LaurentPolynomials
import Jacobians
import MatrixCalculus
import DifferentialFormCalculus

getInput = do
  fh <- openFile "Singular.inp" ReadMode
  conts <- hGetContents fh
  return conts

splitTerms :: String -> [String]
splitTerms ""       = [""]
splitTerms ('+':xs) = []:splitTerms xs
splitTerms ('-':xs) = let ([hd],rst)   = splitAt 1 (splitTerms xs) in "":('-':hd):rst
splitTerms (x:xs)   = let ([hd],rst)   = splitAt 1 (splitTerms xs) in (x:hd):rst

collectNumbers "" = ("","")
collectNumbers (x:xs) = let (nm, rest) = collectNumbers xs in
  if elem x ['0'..'9'] then (x:nm,rest) else ("",x:xs)

colAndConvNum xs = let (str,rst) = collectNumbers xs in
  if str == "" then (1,rst) else (read str :: Int, rst)

parseTerm "" = 1
parseTerm ('-':xs) = (-1) * parseTerm xs
parseTerm xs = Poly [Term coeff [(a,ap)]] * (parseTerm rst) where
  (coeff, rest) = colAndConvNum xs
  terms = if head rest == '*' then tail rest else rest
  a = [head terms]
  (ap,rst) = colAndConvNum $ tail terms


parseLine = sum . map parseTerm . splitTerms

parseInput = do
  inp <- getInput
  return $ map (parseLine) $ lines inp


getJacobian = do
  a <- parseInput
  return $ jacobian (map (\x->[x]) ['a'..'p']) a

sampleComputation = do
  m <- getJacobian
  return $ filter (/=0) . head . matrixAsList $ indMatOnExtAlg m 11
