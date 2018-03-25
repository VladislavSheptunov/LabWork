import unittest
import polynomial as pm

class TestClassPolynomial(unittest.TestCase):

    def setUp(self):
        self.test_srt = str()
        self.tp_A = pm.Polynomial(range(5))
        self.tp_B = pm.Polynomial((1.23, 2.34, 0, -5))
        self.tp_C = pm.Polynomial((0, 0, 0, 0, 0, 0, 2))
        pass

    #@unittest.skip('Temporarily disabled')
    def test_method__eq__(self):
        self.tp_A = pm.Polynomial((1, 2.34, -1, 12))
        self.tp_B = pm.Polynomial((1, 2.34, -1, 12))
        self.assertTrue(self.tp_A == self.tp_B, "Error! Polynomials is not equals")
        self.tp_C = pm.Polynomial((1, 2.3, -1, 23))
        self.assertFalse(self.tp_A == self.tp_C, "Error! Polynomials is equals")

    #@unittest.skip('Temporarily disabled')
    def test_method__str__(self):
        self.test_srt = "10x^10+9x^9+8x^8+7x^7+6x^6+5x^5+4x^4+3x^3+2x^2+x"
        self.tp_A = pm.Polynomial(range(11))
        self.assertEqual(self.tp_A.__str__(), self.test_srt)

        self.test_srt = "-x^2+1.23"
        self.tp_B = pm.Polynomial((1.2345, 0.0001, -1))
        self.assertEqual(self.tp_B.__str__(), self.test_srt)

        self.test_srt = "x^6-x^4+x^2+1"
        self.tp_C = pm.Polynomial((1, 0, 1, 0, -1, 0, 1))
        self.assertEqual(self.tp_C.__str__(), self.test_srt)

        self.test_srt = "-1.1x^2+5.5x-1"
        self.tp_B = pm.Polynomial((1.1, 5.5, -1))
        self.assertNotEqual(self.tp_B.__str__(), self.test_srt)

        self.test_srt = "1x^2+1x-1"
        self.tp_B = pm.Polynomial((-1, 1, 1))
        self.assertNotEqual(self.tp_B.__str__(), self.test_srt)

    #@unittest.skip('Temporarily disabled')
    def test_method__neg__(self):
        self.test_srt = "-4x^6-1.22x^5-44x^3-2.45x+1"
        self.tp_A = pm.Polynomial('-1, 2.45, 0, 44 0, 1.22, 4')
        self.assertEqual((-self.tp_A ).__str__(), self.test_srt)

        self.test_srt = "-x^5-x^3-2x"
        self.tp_A = pm.Polynomial('0, -2, 0, 1 0, -1')
        self.assertNotEqual((-self.tp_A ).__str__(), self.test_srt)

    # @unittest.skip('Temporarily disabled')
    def test_method__pos__(self):
        self.test_srt = "4x^6+1.22x^5+44x^3+2.45x+1"
        self.test_srt = "4x^6+1.22x^5+44x^3+2.45x+1"
        self.tp_A = pm.Polynomial('-1, 2.45, 0, 44 0, -1.22, 4')
        self.assertEqual((+self.tp_A).__str__(), self.test_srt)

        self.test_srt = "-x^5-x^3-2x"
        self.tp_A = pm.Polynomial('0, 2, 0, -1 0, 1')
        self.assertNotEqual((-self.tp_A ).__str__(), self.test_srt)

    #@unittest.skip('Temporarily disabled')
    def test_method__add__(self):
        self.test_srt = "4x^4-2x^3+2x^2+3.34x+1.23"
        self.assertEqual((self.tp_A + self.tp_B).__str__(), self.test_srt)

    #@unittest.skip('Temporarily disabled')
    def test_method__radd__(self):
        self.test_srt = "-5x^3+2.34x+3.33"
        self.assertEqual((2.1 + self.tp_B).__str__(), self.test_srt)

    #@unittest.skip('Temporarily disabled')
    def test_method__iadd__(self):
        self.test_srt = "-5x^3+2x^2+2.34x+2"
        self.tp_A = pm.Polynomial((1.23, 2.34, 0, -5))
        self.tp_A += (0.77, 0, 2)
        self.assertEqual((self.tp_A).__str__(), self.test_srt)

    #@unittest.skip('Temporarily disabled')
    def test_method__sub__(self):
        self.test_srt = "4x^4+8x^3+2x^2-1.34x-1.23"
        self.assertEqual((self.tp_A - self.tp_B).__str__(), self.test_srt)

    #@unittest.skip('Temporarily disabled')
    def test_method__rsub__(self):
        self.test_srt = "5x^3-2.34x+0.87"
        self.assertEqual((2.1 - self.tp_B).__str__(), self.test_srt)

    #@unittest.skip('Temporarily disabled')
    def test_method__isub__(self):
        self.test_srt = "-6x^3+2x^2+2x"
        self.tp_B = pm.Polynomial((1.23, 2.34, 0, -5))
        self.tp_B -= (1.23, 0.34, -2, 1)
        self.assertEqual(self.tp_B.__str__(), self.test_srt)

    #@unittest.skip('Temporarily disabled')
    def test_method__mul__(self):
        self.test_srt = "x^2-1"
        self.tp_A = pm.Polynomial((1, 1))
        self.tp_B = pm.Polynomial((-1, 1))
        self.assertEqual((self.tp_A * self.tp_B).__str__(), self.test_srt)

        self.test_srt = "10x^3-5x^2-10x"
        self.tp_A = pm.Polynomial(-5)
        self.tp_B = pm.Polynomial((0, 2, 1, -2))
        self.assertEqual((self.tp_A * self.tp_B).__str__(),  self.test_srt)

        self.test_srt = "-x^10-2x^9-3x^8"
        self.tp_A = pm.Polynomial((0,0,0,0,0,0,0,0,3,2,1))
        self.tp_B = pm.Polynomial(-1)
        self.assertEqual((self.tp_A * self.tp_B).__str__(), self.test_srt)

        self.tp_B = 0
        self.assertEqual((self.tp_A * self.tp_B).__str__(), '')

        self.test_srt = "x^13-4x^11+5x^5-20x^3"
        self.tp_A = pm.Polynomial((0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 1))
        self.tp_B = pm.Polynomial((0, -4, 0, 1))
        self.assertEqual((self.tp_A * self.tp_B).__str__(), self.test_srt)

    #@unittest.skip('Temporarily disabled')
    def test_method__rmul__(self):
        self.test_srt = "-x^5+3x^3+3x^2-2x-3"
        self.tp_A = pm.Polynomial((3, 2, 0, -1))
        self.tp_B = pm.Polynomial((-1, 0, 1))
        self.assertEqual((self.tp_B * self.tp_A).__str__(), self.test_srt)

    #@unittest.skip('Temporarily disabled')
    def test_method__imul__(self):
        self.test_srt = "-x^5+3x^3+3x^2-2x-3"
        self.tp_A = pm.Polynomial((-1, 0, 1))
        self.tp_A *= (3, 2, 0, -1)
        self.assertEqual((self.tp_A).__str__(), self.test_srt)

    def tearDown(self):
        del self.test_srt
        del self.tp_A
        del self.tp_B
        del self.tp_C
        pass

if __name__ == '__main__':
    unittest.main()
