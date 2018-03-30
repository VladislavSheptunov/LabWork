import re

def RemoveFraction(str_coeff):
    pattern = re.compile('\d\.0\d')
    result = pattern.findall(str_coeff)
    if not result:
        regex = '\.0'
        str = re.sub(regex, '', str_coeff)
    else:
        str = str_coeff
    return str

def IsNumerical(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def ToFloat(_coeffs):
    tmp = []
    for it in _coeffs:
        tmp.append(float(it))
    return tmp

class Polynomial:
    def __init__(self, _coeffs=0):
        self.coeffs = list()
        self.__checkValue(_coeffs)
        self.__addValue(_coeffs)

    def __repr__(self):
       if len(self.coeffs) == 1:
           return 'Polynomial({})'.format(self.coeffs[0])
       else:
           return 'Polynomial({})'.format(tuple(self.coeffs))

    def __str__(self):
        return self.__to_string()

    def __eq__(self, other):
        other = Polynomial(other)
        self.__increase(other)
        zip(self.coeffs, other.coeffs)
        for i, j in zip(self.coeffs, other.coeffs):
            if i != j:
                return False
        return True

    def __ne__(self, other):
        other = Polynomial(other)
        if self.__eq__(other) == False:
            return True
        else:
            return False

    def __neg__(self):
        resPolynomial = Polynomial(self.coeffs)
        for it in range(0, len(resPolynomial.coeffs), 1):
            resPolynomial.coeffs[it] *= -1
        return resPolynomial

    def __pos__(self):
        resPolynomial = Polynomial(self.coeffs)
        for it in range(0, len(resPolynomial.coeffs), 1):
            resPolynomial.coeffs[it] = abs(resPolynomial.coeffs[it])
        return resPolynomial

    def __len__(self):
        return len(self.coeffs)

    def __getitem__(self, key):
        return self.coeffs[key]

    def __iadd__(self, other):
        other = Polynomial(other)
        len_polynomial_result = max(len(self), len(other))
        self.__increase(other)
        for it in range(0, len_polynomial_result, 1):
            self.coeffs[it] += other.coeffs[it]
        return self

    def __add__(self, other):
        other = Polynomial(other)
        resPolynomial = Polynomial(self.coeffs)
        resPolynomial += other
        return resPolynomial

    def __radd__(self, other):
        other = Polynomial(other)
        resPolynomial = Polynomial(self.coeffs)
        other += resPolynomial
        return other

    def __isub__(self, other):
        other = Polynomial(other)
        len_polynomial_result = max(len(self), len(other))
        self.__increase(other)
        for it in range(0, len_polynomial_result, 1):
            self.coeffs[it] -= other.coeffs[it]
        return self

    def __sub__(self, other):
        other = Polynomial(other)
        resPolynomial = Polynomial(self.coeffs)
        resPolynomial -= other
        return resPolynomial

    def __rsub__(self, other):
        other = Polynomial(other)
        resPolynomial = Polynomial(self.coeffs)
        other -= resPolynomial
        return other

    def __imul__(self, other):
        other = Polynomial(other)
        resPolynomial = Polynomial([0 for i in range(len(other) + len(self) - 1)])
        for it1, val1 in enumerate(self.coeffs):
            for it2, val2 in enumerate(other.coeffs):
                resPolynomial.coeffs[it1 + it2] += val1 * val2
        self.coeffs = []
        for it in resPolynomial.coeffs:
            self.coeffs.append(it)
        return self

    def __mul__(self, other):
        other = Polynomial(other)
        resPolynomial = Polynomial(self.coeffs)
        resPolynomial *= other
        return resPolynomial

    def __rmul__(self, other):
        other = Polynomial(other)
        resPolynomial = Polynomial(self.coeffs)
        other *= resPolynomial
        return other

    def __addValue(self, _coeff):
        if (type(_coeff) is int) or (type(_coeff) is float):
            self.coeffs.append(float(_coeff))
        elif type(_coeff) is Polynomial:
            self.coeffs = ToFloat(_coeff.coeffs)
        elif type(_coeff) is str:
            tmp_str = _coeff.split(',')
            for it in range(0, len(tmp_str), 1):
                self.coeffs.append(float(re.sub('\s', '', tmp_str[it])))
        else:
            self.coeffs = ToFloat(list(_coeff))

    def __checkValue(self, _coeff):
        if type(_coeff) is str:
            if len(_coeff) == 1 and IsNumerical(_coeff[0]) == False:
                raise ValueError("Invalid initialization list! Not a number!")
            if (len(_coeff) > 1):
                pattern = re.compile('[^0-9,\s.+-]')
                result = pattern.findall(_coeff)
                if result:
                    raise ValueError("Invalid initialization list! Not a number!")
        if (type(_coeff) is list) or (type(_coeff) is tuple):
            for it in _coeff:
                if (type(it) is not int) and (type(it) is not float):
                    raise ValueError("Error! Invalid type in initialization list")

    def __increase(self, polynomial):
        if len(self) < len(polynomial):
            for it in range(len(self), len(polynomial), 1):
                self.coeffs.insert(0, 0.0)
        else:
            for it in range(len(polynomial), len(self), 1):
                polynomial.coeffs.insert(0, 0.0)

    def __to_string(self, var_string='x', fraction=2):
        res = ''

        if len(self) == self.coeffs.count(0.0):
            return '0'

        first_pow = len(self) - 1
        for i, coeff in enumerate(self.coeffs):
            pow = first_pow - i
            coeff = round(coeff, fraction)

            if coeff:
                if coeff < 0:
                    sign, coeff = '-', -coeff
                elif coeff > 0:
                    sign = '+' if res else ''

                str_coeff = '' if coeff == 1 and pow != 0 else str(coeff)

                if pow == 0:
                    str_power = ''
                elif pow == 1:
                    str_power = var_string
                else:
                    str_power = var_string + '^' + str(pow)

                res += sign + RemoveFraction(str_coeff) + str_power

        return res

    def __del__(self):
        del self.coeffs

if __name__ == "__main__":
    print("For manual testing")
