module SingularPointsTests where

import Control.Monad
import System.IO
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

createAndWriteFiles = do
  forM_ singularFxdPnts (\pnt -> do
    f <- openFile ("SingularFiles/check-"++pnt++".c") WriteMode
    hPutStrLn f (prepend ++ (inputForSingular pnt) ++ append)
    hClose f)

singularFileContent = map ((prepend ++) . (++ append) . inputForSingular) singularFxdPnts

inputForSingular pnt = concat $ zipWith (\n p-> "poly p" ++ show n ++ " = " ++ p ++ ";\n") [1..] polys where
  polys = tail . map f $ localizeNewDefEqs pnt
  f = remOnesAndPluses . remSpace . show . (foldr (.) id ff)
  ff = zipWith (\a b-> plugInIdentity a (Poly [Term 1 [([b],1)]])) (localVars pnt) (['a'..'p'])

remSpace "" = ""
remSpace (' ':xs) = remSpace xs
remSpace (a:xs) = a:remSpace xs

remOnesAndPluses "" = ""
remOnesAndPluses ('1':xs) = remOnesAndPluses xs
remOnesAndPluses ('+':'-':xs) = '-': remOnesAndPluses xs
remOnesAndPluses (a:xs) = a:remOnesAndPluses xs

prepend = "LIB \"primdec.lib\";\noption(redSB);\nring r = 0,(a,b,c,d, e,f,g,h, i,j,k,l, m,n,o,p),dp;\n"
append = "ideal I = p1,p2,p3,p4,p5,p6,p7;\n\nI = groebner(I);\nmatrix H = jacob(I);\nideal J = I,(minor(H,4));\n\nJ = groebner(J);\nJ = radical(J);\n\n\nmatrix M = jacob(J);\nideal K = J,(minor(M,11));\nK = groebner (K);\nprint(J);\nprint(K);\n"
