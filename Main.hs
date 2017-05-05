
import MatrixCalculus
import Jacobians
import EigenvalueList

main = mapM_ (putStrLn) . 
  zipWith (\v m -> "\\wil J_{"++ v ++ "} = \n" ++ m) eigVecsInXmin $ 
    map (latexPrintMatrix . jacobianOfNewDefEqsAtO) eigVecsInXmin
