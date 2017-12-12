
import System.IO
import MatrixCalculus
import LieG2Basis

m = map (latexPrintMatrix . (twiceIdentity *)) g2TildeBasis

m' = map (\(n,mt) -> "$$\n2\\tilde X_{"++show n++ "}="++ mt ++"$$\n\n") $ zip [1..] m

g = map latexPrintMatrix g2Basis

g' = map (\(n,mt) -> "$$\nX_{"++show n++ "}="++ mt ++"$$\n\n") $ zip [1..] g

-- g' is basis for g2 algebra in "standard" basis
-- m' is basis for g2 algebra in transformed basis
--    first, everything in g' is transformed and
--    then multiplied by 2
main = do
  putStrLn $ latexPrintMatrix singleMatParabolicStandardBasis
--  mapM_ putStrLn g'
