from numpy.polynomial import Polynomial
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

def FormMonom(coeff, frac):
    if coeff < 0:
        if abs(coeff) == 1:
            str_tmp = ' - ' + 'x'
        else:
            str_tmp = ' - ' + str(abs(round(coeff, frac))) + 'x'
    elif coeff > 0:
        if abs(coeff) == 1:
            str_tmp = ' + ' + 'x'
        else:
            str_tmp = ' + ' + str(abs(round(coeff, frac))) + 'x'
    else:
        str_tmp = "empty"
    return str_tmp

def GetMonom(monom_number, coeff, code, frac):
    if monom_number == 0:
        if coeff < 0:
            str_tmp = ' - ' + str(abs(round(coeff, frac)))
        elif coeff > 0:
            str_tmp = ' + ' + str(abs(round(coeff, frac)))
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
    arg.pop(0)
    if arg[0] == '+':
        arg.pop(0)
        arg.pop(0)
    else:
        arg.pop(1)
    for it in arg:
        ret_str += str(it)
    return ret_str

def GetPolinomial(list_monom):
    ret_str = str()
    for it in range(0, len(list_monom), 1):
        ret_str += str(list_monom[it])
    return ret_str

def IsFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def ToFloat(coeff):
    tmp = []
    for it in coeff:
        tmp.append(float(it))
    return tmp

class Polynomial:
    code = {'0': chr(8304), '1': chr(185), '2': chr(178),
            '3': chr(179), '4': chr(8308), '5': chr(8309),
            '6': chr(8310), '7': chr(8311), '8': chr(8312),
            '9': chr(8313)}

    def __init__(self, coeff = 0):
        self._coeff = list()
        self.__checkValue(coeff)
        self.__addValue(coeff)

    def __repr__(self):
       if len(self._coeff) == 1:
            return 'Polynomial({})'.format(self._coeff[0])
       else:
           tmp_list = self._coeff
           for it in reversed(tmp_list):
               if it == 0:
                   tmp_list.remove(it)
               else:
                   break
           return 'Polynomial({})'.format(tuple(tmp_list))

    def __str__(self):
        list_monom = list()
        len_poly = len(self._coeff)
        if len_poly > 1 and len_poly - self._coeff.count(0) != 0:
            for monom_number in range(0, len_poly, 1):
                tmp_str = GetMonom(monom_number, self._coeff[monom_number], self.code, 3)
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
        if len(self._coeff) != len(other._coeff):
            return False
        zip(self._coeff, other._coeff)
        for i, j in zip(self._coeff, other._coeff):
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
        resPolynomial = Polynomial(self._coeff)
        for it in range(0, len(resPolynomial._coeff),1):
            resPolynomial._coeff[it] *= -1
        return resPolynomial

    def __pos__(self):
        resPolynomial = Polynomial(self._coeff)
        for it in range(0, len(resPolynomial._coeff),1):
            resPolynomial._coeff[it] = abs(resPolynomial._coeff[it])
        return resPolynomial

    def __bool__(self):
        if len(self._coeff) == 1 and self._coeff[0] == 0:
            return False
        else:
            return True

    def __iadd__(self, other):
        other = Polynomial(other)
        len_poly = len(other._coeff)
        if len(self._coeff) <= len_poly:
            for it in range(len(self._coeff),len(self._coeff) + (len_poly - len(self._coeff)), 1):
                self._coeff.append(0)
        for it in range(0,len_poly,1):
            self._coeff[it] += other._coeff[it]
            tmp = self._coeff[it]
            if tmp - self._coeff[it] == 0:
                int(self._coeff[it])
        return self

    def __add__(self, other):
        other = Polynomial(other)
        resPolynomial = Polynomial(self._coeff)
        resPolynomial += other
        return resPolynomial

    def __radd__(self, other):
        other = Polynomial(other)
        resPolynomial = Polynomial(self._coeff)
        other += resPolynomial
        return other

    def __isub__(self, other):
        other = Polynomial(other)
        len_poly = len(other._coeff)
        if len(self._coeff) <= len_poly:
            for it in range(len(self._coeff),len(self._coeff) + (len_poly - len(self._coeff)), 1):
                self._coeff.append(0)
        for it in range(0,len_poly,1):
            self._coeff[it] -= other._coeff[it]
        return self

    def __sub__(self, other):
        other = Polynomial(other)
        resPolynomial = Polynomial(self._coeff)
        resPolynomial -= other
        return resPolynomial

    def __rsub__(self, other):
        other = Polynomial(other)
        resPolynomial = Polynomial(self._coeff)
        other -= resPolynomial
        return other

    def __imul__(self, other):
        other = Polynomial(other)
        len_poly = len(self._coeff) + len(other._coeff)
        resPolynomial = Polynomial([0 for i in range(len_poly)])
        for i in range(0, len(self._coeff), 1):
            for j in range(0, len(other._coeff), 1):
                resPolynomial._coeff[i + j] += self._coeff[i] * other._coeff[j]
        self._coeff = []
        for it in resPolynomial._coeff:
            self._coeff.append(it)
        return self


    def __mul__(self, other):
        other = Polynomial(other)
        resPolynomial = Polynomial(self._coeff)
        resPolynomial *= other
        return resPolynomial

    def __rmul__(self, other):
        other = Polynomial(other)
        resPolynomial = Polynomial(self._coeff)
        other *= resPolynomial
        return other

    def __addValue(self, coeff):
        if (type(coeff) is int) or (type(coeff) is float):
            self._coeff.append(float(coeff))
        elif type(coeff) is Polynomial:
            self._coeff = ToFloat(coeff._coeff)
        elif type(coeff) is str:
            tmp_str = re.sub('[\,]', '', coeff).split(' ')
            for it in range(0,len(tmp_str),1):
                if IsFloat(tmp_str[it]):
                    self._coeff.append(float(tmp_str[it]))
        else:
            self._coeff = ToFloat(list(coeff))

    def __checkValue(self, coeff):
        if type(coeff) is str:
            if len(coeff) == 1 and IsFloat(coeff[0]) == False:
                raise ValueError("Invalid initialization list! Not a number!")
            if (len(coeff) > 1):
                tmp_checkList = re.sub('[\,]', '', coeff).split(' ')
                for it in tmp_checkList:
                    if IsFloat(it) == False:
                        raise ValueError("Invalid initialization list! Not a number!")
        if (type(coeff) is list) or (type(coeff) is tuple):
            for it in coeff:
                if (type(it) is not int) and (type(it) is not float):
                    raise ValueError("Error! Invalid type in initialization list")

    def __del__(self):
        del self._coeff


if __name__ == "__main__":
    print("For manual testing")
