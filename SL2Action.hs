module SL2Actions where

import ComplexNumbers (i)
import ComplexOctonions
import LaurentPolynomials

t1 = Term  1 [("a",1)]
t2 = Term oi [("b",1)]
t3 = Term oj [("c",1)]
t4 = Term ok [("d",1)]

g = Poly [t1, t2, t3, t4]

ix = scalarMult i
te0 = 1 + ix oi
te1 = 1 - ix oi
te2 = oj + ix ok
te3 = oj - ix ok
te4 = ol + ix oli
te5 = ol - ix oli
te6 = olj + ix olk
te7 = olj - ix olk
cte0 = conv te0
cte1 = conv te1
cte2 = conv te2
cte3 = conv te3
cte4 = conv te4
cte5 = conv te5
cte6 = conv te6
cte7 = conv te7

tBasis = [te0, te1, te2, te3,
          te4, te5, te6, te7]

convertedtBasis = map conv tBasis

conv :: Octonion -> LPoly Octonion
conv o = Poly [Term o []]

ml = -ol
c1 = conv 1
ci = conv oi
cj = conv oj
ck = conv ok
cl = conv ol
cli = conv oli
clj = conv olj
clk = conv olk

