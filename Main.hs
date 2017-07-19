
import System.IO
import MatrixCalculus
import Jacobians
import SingularTest

theMat = indMatOnExtAlg (rows [1..7] . jacobianInLocalCoords "0123" $ localizeNewDefEqs "0123") 4

main = do 
  hSetBuffering stdout NoBuffering
  sc <- sampleComputation
  putStrLn $ show sc
