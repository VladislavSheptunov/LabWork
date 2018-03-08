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
        self.test_srt = "{}x² {} {}x {} 1".format('- ', '+', '2.34', '+')
        P = pm.Polynomial((1,2.34,-1,0))
        self.assertEqual(P.__str__(), self.test_srt)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
