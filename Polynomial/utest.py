import unittest
import polynomial as pm

class TestClassPolynomial(unittest.TestCase):
    def setUp(self):
        self.tp_A = pm.Polynomial((1, 2, 0, -3))
        self.tp_B = pm.Polynomial((1.23, 2.34, 0, -5))
        self.tp_C = pm.Polynomial((0, 0, 0, 0, 0, 0, 2))
        self.tp_tmp = pm.Polynomial()
        pass

    #@unittest.skip('Temporarily disabled')
    def test_method__eq__(self):
        self.tp_A = pm.Polynomial((1, 2.34, -1, 12))
        self.tp_B = pm.Polynomial((1, 2.34, -1, 12))
        self.assertTrue(self.tp_A == self.tp_B, "Error! Polynomials is not equals")
        self.tp_C = pm.Polynomial((1, 2.3, -1, 23))
        self.assertFalse(self.tp_A == self.tp_C, "Error! Polynomials is equals")

        self.tp_tmp = pm.Polynomial([0, 0, 0])
        self.assertTrue(self.tp_tmp == pm.Polynomial([0]), "Error! 0 == 0")

        self.tp_tmp = pm.Polynomial([0, 0, 3, 5])
        self.assertTrue(self.tp_tmp == pm.Polynomial([3, 5]), "Error! 0x^3 + 0x^2 + 3x + 5 == 3x + 5")

    #@unittest.skip('Temporarily disabled')
    def test_method__str__positive(self):
        self.tp_A = pm.Polynomial(range(11))
        self.assertEqual(self.tp_A.__str__(), "x^9+2x^8+3x^7+4x^6+5x^5+6x^4+7x^3+8x^2+9x+10")

        self.tp_B = pm.Polynomial((1.2345, 0.0001, -1))
        self.assertEqual(self.tp_B.__str__(), "1.23x^2-1")

        self.tp_C = pm.Polynomial((1, 0, 1, 0, -1, 0, 1))
        self.assertEqual(self.tp_C.__str__(), "x^6+x^4-x^2+1")

        self.tp_A = pm.Polynomial((1.2, 1.0001, -1))
        self.assertEqual(self.tp_A.__str__(), "1.2x^2+x-1")

        self.tp_B = pm.Polynomial((1.2, 1.0101, -1))
        self.assertEqual(self.tp_B.__str__(), "1.2x^2+1.01x-1")

        self.tp_C = pm.Polynomial()
        self.assertEqual(self.tp_C.__str__(), "0")

        self.tp_A = pm.Polynomial(3)
        self.assertEqual(self.tp_A.__str__(), "3")

    #@unittest.skip('Temporarily disabled')
    def test_method__str__negative(self):
        self.tp_A = pm.Polynomial((1, 5.5, -1))
        self.assertNotEqual(self.tp_A.__str__(), "1x^2+5.5x-1")

        self.tp_B = pm.Polynomial((1, 5.0, -1))
        self.assertNotEqual(self.tp_A.__str__(), "x^2+5.0x-1")

        self.tp_C = pm.Polynomial((1, 5.001, 1))
        self.assertNotEqual(self.tp_A.__str__(), "x^2+5.001x+1")

    #@unittest.skip('Temporarily disabled')
    def test_method__neg__positive(self):
        self.assertEqual((-self.tp_B).__str__(), "-1.23x^3-2.34x^2+5")

    #@unittest.skip('Temporarily disabled')
    def test_method__neg__negative(self):
        self.assertNotEqual((-self.tp_B).__str__(), "-1.23x^3-2.34x^2-5")
        self.assertNotEqual((-self.tp_B).__str__(), "1.23x^3+2.34x^2+5")

    #@unittest.skip('Temporarily disabled')
    def test_method__pos__positive(self):
        self.tp_A = pm.Polynomial('-1, 2.45, 0, 44, 0, -1.22, 4')
        self.tp_B = pm.Polynomial('-1, -2, 3, -4 , 5')
        self.assertEqual((+self.tp_A).__str__(), "x^6+2.45x^5+44x^3+1.22x+4")
        self.assertEqual((+self.tp_B).__str__(), "x^4+2x^3+3x^2+4x+5")

    #@unittest.skip('Temporarily disabled')
    def test_method__pos__negative(self):
        self.assertNotEqual((-self.tp_B).__str__(), "1.23x^3+2.34x^2-5")
        self.assertNotEqual((-self.tp_B).__str__(), "+1.23x^3+2.34x^2+5")

    #@unittest.skip('Temporarily disabled')
    def test_method__iadd__positive(self):
        self.tp_A += (0.77, 1, 2)
        self.assertEqual((self.tp_A).__str__(), "x^3+2.77x^2+x-1")

        self.tp_A += (1, 3, 0.23, 2, 1)
        self.assertEqual((self.tp_A).__str__(), "x^4+4x^3+3x^2+3x")

        self.tp_A += pm.Polynomial("-3, 1")
        self.assertEqual((self.tp_A).__str__(), "x^4+4x^3+3x^2+1")

        self.tp_tmp += 1
        self.assertEqual((self.tp_tmp).__str__(), "1")

    #@unittest.skip('Temporarily disabled')
    def test_method__iadd__negative(self):
        self.tp_A += pm.Polynomial([-2, 0, 3])
        self.assertNotEqual((self.tp_A).__str__(), "1")

        self.tp_B += (0.77, -2.34, 0, 5)
        self.assertNotEqual((self.tp_B).__str__(), "2.0x^4")

        self.tp_tmp += 1.0
        self.assertNotEqual((self.tp_tmp).__str__(), "1.0")

    #@unittest.skip('Temporarily disabled')
    def test_method__add__(self):
        self.assertEqual((self.tp_A + self.tp_B).__str__(), "2.23x^3+4.34x^2-8")

        self.assertEqual((self.tp_A + 3).__str__(), "x^3+2x^2")

        self.tp_tmp = pm.Polynomial("-1, -2, 0, 4")
        self.assertEqual((self.tp_tmp + self.tp_A).__str__(), "1")

    #@unittest.skip('Temporarily disabled')
    def test_method__radd__(self):
        self.assertEqual((2.1 + self.tp_A).__str__(), "x^3+2x^2-0.9")

        self.assertEqual((self.tp_B + self.tp_A).__str__(), "2.23x^3+4.34x^2-8")

        self.assertEqual((3 + self.tp_A).__str__(), "x^3+2x^2")

        self.tp_tmp = pm.Polynomial("-1, -2, 0, 3")
        self.assertEqual((self.tp_A + self.tp_tmp).__str__(), "0")

    #@unittest.skip('Temporarily disabled')
    def test_method__isub__positive(self):
        self.tp_B -= (1.23, 0.34, -2, 1)
        self.assertEqual(self.tp_B.__str__(), "2x^2+2x-6")

        self.tp_tmp = pm.Polynomial(1)
        self.tp_tmp -= 1
        self.assertEqual(self.tp_tmp.__str__(), "0")

    #@unittest.skip('Temporarily disabled')
    def test_method__isub__negative(self):
        self.tp_B -= (1.23, 0.34, -2, 1)
        self.assertNotEqual(self.tp_B.__str__(), "-6x^3+2x^2+2x")

        self.tp_tmp -= 1.0
        self.assertNotEqual((self.tp_tmp).__str__(), "1")

    #@unittest.skip('Temporarily disabled')
    def test_method__sub__(self):
        self.tp_A = pm.Polynomial("3, 3, 3")
        self.tp_B = pm.Polynomial("2, 2, 2")
        self.assertEqual((self.tp_A - self.tp_B).__str__(), "x^2+x+1")

        self.tp_C = pm.Polynomial("3, 3")
        self.assertEqual((self.tp_A - self.tp_C).__str__(), "3x^2")

        self.assertEqual((self.tp_A - [-1, 3, 3, 3]).__str__(), "x^3")

    #@unittest.skip('Temporarily disabled')
    def test_method__rsub__(self):
        self.tp_A = pm.Polynomial((1, 2, 0, -3))
        self.assertEqual((2.1 - self.tp_A).__str__(), "-x^3-2x^2+5.1")

    #@unittest.skip('Temporarily disabled')
    def test_method__imul__(self):
        self.tp_A = pm.Polynomial((1, 0, -1))
        self.tp_A *= (-1, 0, 2, 3)
        self.assertEqual((self.tp_A).__str__(), "-x^5+3x^3+3x^2-2x-3")
        self.tp_A *= -1
        self.assertEqual((self.tp_A).__str__(), "x^5-3x^3-3x^2+2x+3")

    #@unittest.skip('Temporarily disabled')
    def test_method__mul__(self):
        self.tp_A = pm.Polynomial((1, 1))
        self.tp_B = pm.Polynomial((1, -1))
        self.assertEqual((self.tp_A * self.tp_B).__str__(), "x^2-1")

        self.tp_A = pm.Polynomial(-5)
        self.tp_B = pm.Polynomial((-2, 1, 2, 0))
        self.assertEqual((self.tp_A * self.tp_B).__str__(), "10x^3-5x^2-10x")

        self.tp_A = pm.Polynomial((1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0))
        self.tp_B = pm.Polynomial(-1)
        self.assertEqual((self.tp_A * self.tp_B).__str__(), "-x^10-2x^9-3x^8")

        self.tp_B = 0
        self.assertEqual((self.tp_A * self.tp_B).__str__(), '0')

        self.tp_A = pm.Polynomial((1, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0))
        self.tp_B = pm.Polynomial((1, 0, -4, 0))
        self.assertEqual((self.tp_A * self.tp_B).__str__(), "x^13-4x^11+5x^5-20x^3")

    #@unittest.skip('Temporarily disabled')
    def test_method__rmul__(self):
        self.tp_A = pm.Polynomial((-1, 0, 2, 3))
        self.tp_B = pm.Polynomial((1, 0, -1))
        self.assertEqual((self.tp_B * self.tp_A).__str__(), "-x^5+3x^3+3x^2-2x-3")

        self.assertEqual((-1 * self.tp_B).__str__(), "-x^2+1")

    def tearDown(self):
        del self.tp_A
        del self.tp_B
        del self.tp_C
        del self.tp_tmp
        pass

if __name__ == '__main__':
    unittest.main()
