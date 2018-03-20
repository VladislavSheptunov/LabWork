import re

def RemoveFraction(monom, number_monom):
    if number_monom == 0:
        pattern = re.compile('\d\.0\d')
        result = pattern.findall(monom)
        if not result:
            regex = '\.0'
            str = re.sub(regex, '', monom)
        else:
            str = monom
    else:
        pattern = re.compile('\d.0x')
        result = pattern.findall(monom)
        if result:
            regex = '\.0'
            str = re.sub(regex, '', monom)
        else:
            str = monom
    return str

def SetItem(_coeffs, frac, sig):
    if abs(_coeffs) == 1:
        return sig + 'x'
    else:
        return sig + str(abs(round(_coeffs, frac))) + 'x'

def FormMonom(_coeffs, frac):
    if _coeffs < 0:
        return SetItem(_coeffs, frac, '-')
    elif _coeffs > 0:
        return SetItem(_coeffs, frac, '+')
    else:
        str_tmp = "empty"
    return str_tmp

def GetMonom(monom_number, coeff, code, frac):
    if monom_number == 0:
        if coeff < 0:
            str_tmp = '-' + str(abs(round(coeff, frac)))
        elif coeff > 0:
            str_tmp = '+' + str(abs(round(coeff, frac)))
        else:
            str_tmp = "empty"
    if monom_number == 1:
        str_tmp = FormMonom(coeff, frac)
    if monom_number > 1:
        str_tmp = FormMonom(coeff, frac)
        if str_tmp == "empty":
            return "empty"
        tmp = map(int, str(monom_number))
        for jt in tmp:
            str_tmp += str(code[str(jt)])
        del tmp
    if str_tmp != "empty":
        str_tmp = RemoveFraction(str_tmp, monom_number)
    return str_tmp

def EditSeniorMonom(arg):
    ret_str = str()
    if arg[0] == '+':
        arg.pop(0)
    for it in arg:
        ret_str += str(it)
    return ret_str

def GetPolinomial(list_monom):
    ret_str = str()
    for it in range(0, len(list_monom), 1):
        ret_str += str(list_monom[it])
    return ret_str

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
    code = {'0': chr(8304), '1': chr(185), '2': chr(178),
            '3': chr(179), '4': chr(8308), '5': chr(8309),
            '6': chr(8310), '7': chr(8311), '8': chr(8312),
            '9': chr(8313)}

    def __init__(self, _coeffs = 0):
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
        list_monom = list()
        len_poly = len(self.coeffs)
        if len_poly > 1 and len_poly - self.coeffs.count(0) != 0:
            for monom_number in range(0, len_poly, 1):
                tmp_str = GetMonom(monom_number, self.coeffs[monom_number], self.code, 3)
                if tmp_str != "empty":
                    list_monom.append(tmp_str)
                else:
                    continue
            list_monom.reverse()
            list_monom[0] = EditSeniorMonom(list(list_monom[0]))
        else:
            list_monom.append(0)
        return GetPolinomial(list_monom)

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

    def __bool__(self):
        if len(self.coeffs) == 1 and self.coeffs[0] == 0:
            return False
        else:
            return True

    def __iadd__(self, other):
        other = Polynomial(other)
        len_poly = len(other.coeffs)
        if len(self.coeffs) <= len_poly:
            for it in range(len(self.coeffs), len(self.coeffs) + (len_poly - len(self.coeffs)), 1):
                self.coeffs.append(0)
        for it in range(0,len_poly,1):
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
        for it in range(0,len_poly,1):
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


if __name__ == "__main__":
    print("For manual testing")
