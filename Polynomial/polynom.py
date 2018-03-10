from numpy.polynomial import Polynomial
import re

def GetMonom(it, _coeff, code):
    str_tmp = str()
    if it == 0:
        if _coeff[it] < 0:
            str_tmp = ' - ' + str(abs(round(_coeff[it], 2)))
        elif _coeff[it] > 0:
            str_tmp = ' + ' + str(abs(round(_coeff[it], 2)))
        else:
            return "continue"
    if it == 1:
        if _coeff[it] < 0:
            if abs(_coeff[it]) == 1:
                str_tmp = ' - ' + 'x'
            else:
                str_tmp = ' - ' + str(abs(round(_coeff[it], 2))) + 'x'
        elif _coeff[it] > 0:
            if abs(_coeff[it]) == 1:
                str_tmp = ' + ' + 'x'
            else:
                str_tmp = ' + ' + str(abs(round(_coeff[it], 2))) + 'x'
        else:
            return "continue"
    if it > 1:
        if _coeff[it] < 0:
            if abs(_coeff[it]) == 1:
                str_tmp += ' - ' + 'x'
            else:
                str_tmp += ' - ' + str(abs(round(_coeff[it], 2))) + 'x'
        elif _coeff[it] > 0:
            if abs(_coeff[it]) == 1:
                str_tmp += ' + ' + 'x'
            else:
                str_tmp += ' + ' + str(abs(round(_coeff[it], 2))) + 'x'
        else:
            return "continue"
        tmp = map(int, str(it))
        for jt in tmp:
            str_tmp += str(code[str(jt)])
        del tmp
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

def isint(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def isfloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

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
            for it in range(0, len_poly, 1):
                tmp_str = GetMonom(it, self._coeff, self.code)
                if tmp_str == "continue":
                    continue
                else:
                    list_monom.append(tmp_str)
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

    def __addValue(self, coeff):
        if (type(coeff) is int) or (type(coeff) is float):
            self._coeff.append(coeff)
        elif type(coeff) is Polynomial:
            self._coeff = coeff._coeff
        elif type(coeff) is str:
            tmp_str = re.sub('[\,]', '', coeff).split(' ')
            for it in range(0,len(tmp_str),1):
                if isint(tmp_str[it]):
                    self._coeff.append(int(tmp_str[it]))
                elif isfloat(tmp_str[it]):
                    self._coeff.append(float(tmp_str[it]))
        else:
            self._coeff = list(coeff)

    def __checkValue(self, coeff):
        if type(coeff) is str:
            if len(coeff) == 1 and isfloat(coeff[0]) == False:
                raise ValueError("Invalid initialization list! Not a number!")
            if (len(coeff) > 1):
                tmp_checkList = re.sub('[\,]', '', coeff).split(' ')
                for it in tmp_checkList:
                    if isfloat(it) == False:
                        raise ValueError("Invalid initialization list! Not a number!")
        if (type(coeff) is list) or (type(coeff) is tuple):
            for it in coeff:
                if (type(it) is not int) and (type(it) is not float):
                    raise ValueError("Error! Invalid type in initialization list")

    def __del__(self):
        del self._coeff

# ========================================================================================================#
if __name__ == "__main__":
    t_P = Polynomial((1.23, 2.34, 0, -5))
    t_P += (0, 77, 0, 2)
    print(t_P)
