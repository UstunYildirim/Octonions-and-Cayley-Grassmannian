
import System.IO
import MatrixCalculus
import LieG2Basis

m = map (latexPrintMatrix . (twiceIdentity *)) g2TildeBasis

m' = map (\(n,mt) -> "$$\nX_{"++show n++ "}="++ mt ++"$$\n\n") $ zip [1..] m

main = do
  mapM_ putStrLn m'
