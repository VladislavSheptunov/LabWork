import re

class Polynomial:
    def __init__(self, _coeffs=0):
        self.coeffs = list()
        self.__checkValue(_coeffs)
        self.__addValue(_coeffs)

    def __repr__(self):
       if len(self.coeffs) == 1:
            return 'Polynomial({})'.format(self.coeffs[0])
       else:
           tmp_list = self.coeffs
           for it in reversed(tmp_list):
               if it == 0:
                   tmp_list.remove(it)
               else:
                   break
           return 'Polynomial({})'.format(tuple(tmp_list))

    def __str__(self):
        polynom = self.coeffs
        polynom.reverse()
        return FormStringPolynomial(polynom)

    def __eq__(self, other):
        other = Polynomial(other)
        if len(self.coeffs) != len(other.coeffs):
            return False
        zip(self.coeffs, other.coeffs)
        for i, j in zip(self.coeffs, other.coeffs):
            if i != j:
                return False
        return True

    def __ne__(self, other):
        other = Polynomial(other)
        if self.__eq__(other) == False:
            return True;
        else:
            return False;

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
        len_poly = len(other.coeffs)
        if len(self.coeffs) <= len_poly:
            for it in range(len(self.coeffs), len(self.coeffs) + (len_poly - len(self.coeffs)), 1):
                self.coeffs.append(0)
        for it in range(0, len_poly, 1):
            self.coeffs[it] += other.coeffs[it]
            tmp = self.coeffs[it]
            if tmp - self.coeffs[it] == 0:
                int(self.coeffs[it])
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
        len_poly = len(other.coeffs)
        if len(self.coeffs) <= len_poly:
            for it in range(len(self.coeffs), len(self.coeffs) + (len_poly - len(self.coeffs)), 1):
                self.coeffs.append(0)
        for it in range(0, len_poly, 1):
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
        len_poly = len(self.coeffs) + len(other.coeffs)
        resPolynomial = Polynomial([0 for i in range(len_poly)])
        for i in range(0, len(self.coeffs), 1):
            for j in range(0, len(other.coeffs), 1):
                resPolynomial.coeffs[i + j] += self.coeffs[i] * other.coeffs[j]
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
            tmp_str = re.sub('[\,]', '', _coeff).split(' ')
            for it in range(0,len(tmp_str),1):
                if IsNumerical(tmp_str[it]):
                    self.coeffs.append(float(tmp_str[it]))
        else:
            self.coeffs = ToFloat(list(_coeff))

    def __checkValue(self, _coeff):
        if type(_coeff) is str:
            if len(_coeff) == 1 and IsNumerical(_coeff[0]) == False:
                raise ValueError("Invalid initialization list! Not a number!")
            if (len(_coeff) > 1):
                tmp_checkList = re.sub('[\,]', '', _coeff).split(' ')
                for it in tmp_checkList:
                    if IsNumerical(it) == False:
                        raise ValueError("Invalid initialization list! Not a number!")
        if (type(_coeff) is list) or (type(_coeff) is tuple):
            for it in _coeff:
                if (type(it) is not int) and (type(it) is not float):
                    raise ValueError("Error! Invalid type in initialization list")

    def __del__(self):
        del self.coeffs

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

def FormStringPolynomial(polynomial, var_string='x'):
    res = ''
    first_pow = len(polynomial) - 1
    for i, coeff in enumerate(polynomial):
        power = first_pow - i
        coeff = round(coeff, 2)

        if coeff:
            if coeff < 0:
                sign, coeff = '-', -coeff
            elif coeff > 0:
                sign = ('+' if res else '')

            str_coeff = '' if coeff == 1 and power != 0 else str(coeff)

            if power == 0:
                str_power = ''
            elif power == 1:
                str_power = var_string
            else:
                str_power = var_string + '^' + str(power)

            res += sign + RemoveFraction(str_coeff) + str_power
    return res


def ToFloat(_coeffs):
    tmp = []
    for it in _coeffs:
        tmp.append(float(it))
    return tmp

if __name__ == "__main__":
    print("For manual testing")
    #print(Polynomial("1, 2, 3"))
    #print(Polynomial("1.001, 2, 3"))
