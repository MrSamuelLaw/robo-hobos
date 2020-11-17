import unittest
import sympy as sym
from mathtools import symbolic


class test_symbolic(unittest.TestCase):

    def test_sym2poly(self):
        # setup the test variables
        s, K = sym.symbols("s, K")
        oltf = 1/(s+1)
        cltf = oltf/(1 + oltf)
        num, den = sym.fraction(cltf)

        # num = 1, test that it can handle
        # an expression with no symbols
        res1 = symbolic.sym2poly(num)
        res2 = symbolic.sym2poly(num, var=s)
        self.assertEqual(res1, res2)

        # den = s+2, test that it can handle a
        # mono variable expression
        res1 = symbolic.sym2poly(den)
        res2 = symbolic.sym2poly(den, var=s)
        self.assertEqual(res1, res2)

        # assert that you have to specify variable
        # if using symbolic coefficients
        p = K*s**2 + 5*s
        with self.assertRaises(Exception):
            symbolic.sym2poly(p)
        symbolic.sym2poly(p, var=s)

        # asser that numerator and denominator
        # can be processed at the same time
        numden = symbolic.sym2poly(num, den)
        self.assertEqual(len(numden), 2)

    def test_RHtable(self):
        pass


if __name__ == "__main__":
    unittest.main()