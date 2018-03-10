import unittest
import polynom as pm


class TestClassPolynomial(unittest.TestCase):

    def setUp(self):
        pass

    def test_method__eq__(self):
        P1 = pm.Polynomial((1,2.34,-1,12))
        P2 = pm.Polynomial((1, 2.34, -1, 12))
        self.assertTrue(P1 == P2, "Error! Polynomials is not equals")
        P3 = pm.Polynomial((1, 2.3, -1, 23))
        self.assertFalse(P1 == P3, "Error! Polynomials is equals")

    def test_method__str__(self):
        self.test_srt1 = "10x¹⁰ + 9x⁹ + 8x⁸ + 7x⁷ + 6x⁶ + 5x⁵ + 4x⁴ + 3x³ + 2x² + x"
        self.test_srt2 = "-x² + 0.0x + 1.23"
        self.test_srt3 = "x⁶ - x⁴ + x² + 1"
        t_P1 = pm.Polynomial(range(11))
        t_P2 = pm.Polynomial((1.2345, 0.0001, -1))
        t_P3 = pm.Polynomial((1, 0, 1, 0, -1, 0, 1))

        self.assertEqual(t_P1.__str__(), self.test_srt1)
        self.assertEqual(t_P2.__str__(), self.test_srt2)
        self.assertEqual(t_P3.__str__(), self.test_srt3)

    def test_method__neg__(self):
        self.test_srt = "-4x⁶ - 1.22x⁵ - 44x³ - 2.45x + 1"
        t_P = pm.Polynomial('-1, 2.45, 0, 44 0, 1.22, 4')
        self.assertEqual((-t_P).__str__(), self.test_srt)

    def test_method__pos__(self):
        self.test_srt = "4x⁶ + 1.22x⁵ + 44x³ + 2.45x + 1"
        t_P = pm.Polynomial('-1, 2.45, 0, 44 0, -1.22, 4')
        self.assertEqual((+t_P).__str__(), self.test_srt)

    def test_method__bool__(self):
        t_P1 = pm.Polynomial()
        self.assertFalse(bool(t_P1))
        t_P2 = pm.Polynomial((0, 4.3, 9.2))
        self.assertTrue(bool(t_P2))

    def test_method__add__(self):
        self.test_srt = "4x⁴ - 2x³ + 2x² + 3.34x + 1.23"
        t_P1 = pm.Polynomial(range(5))
        t_P2 = pm.Polynomial((1.23, 2.34, 0, -5))
        self.assertEqual((t_P1 + t_P2).__str__(), self.test_srt)

    def test_method__radd__(self):
        self.test_srt = "-5x³ + 2.34x + 3.33"
        t_P = pm.Polynomial((1.23, 2.34, 0, -5))
        self.assertEqual((2.1 + t_P).__str__(), self.test_srt)

    def test_method__iadd__(self):
        self.test_srt = "-5x³ + 2x² + 2.34x + 2.0"
        t_P = pm.Polynomial((1.23, 2.34, 0, -5))
        t_P += (0.77, 0, 2)
        self.assertEqual((t_P).__str__(), self.test_srt)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
