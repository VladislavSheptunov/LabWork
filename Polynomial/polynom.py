from numpy.polynomial import Polynomial

def GetMonom(it, _coeff, code):
    str_tmp = str()
    if it == 0:
        if _coeff[it] < 0:
            str_tmp = ' - ' + str(abs(_coeff[it]))
        elif _coeff[it] > 0:
            str_tmp = ' + ' + str(abs(_coeff[it]))
        else:
            return "continue"
    if it == 1:
        if _coeff[it] < 0:
            if abs(_coeff[it]) == 1:
                str_tmp = ' - ' + 'x'
            else:
                str_tmp = ' - ' + str(abs(_coeff[it])) + 'x'
        elif _coeff[it] > 0:
            if abs(_coeff[it]) == 1:
                str_tmp = ' + ' + 'x'
            else:
                str_tmp = ' + ' + str(abs(_coeff[it])) + 'x'
        else:
            return "continue"
    if it > 1:
        if _coeff[it] < 0:
            if abs(_coeff[it]) == 1:
                str_tmp += ' - ' + 'x'
            else:
                str_tmp += ' - ' + str(abs(_coeff[it])) + 'x'
        elif _coeff[it] > 0:
            if abs(_coeff[it]) == 1:
                str_tmp += ' + ' + 'x'
            else:
                str_tmp += ' + ' + str(abs(_coeff[it])) + 'x'
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
    for it in arg:
        ret_str += str(it)
    return ret_str

def GetPolinomial(list_monom):
    ret_str = str()
    for it in range(0, len(list_monom), 1):
        ret_str += str(list_monom[it])
    return ret_str


class Polynomial:
    code = {'0': chr(8304), '1': chr(185), '2': chr(178),
            '3': chr(179), '4': chr(8308), '5': chr(8309),
            '6': chr(8310), '7': chr(8311), '8': chr(8312),
            '9': chr(8313)}

    def __init__(self, coeff = 0):
        self._coeff = list()
        self.__checkValue(coeff)
        self.__addValue(coeff)
        self._coeff = tuple(self._coeff)

    def __repr__(self):
       if len(self._coeff) == 1:
            return 'Polynom({})'.format(self._coeff[0])
       else:
            return 'Polynom({})'.format(self._coeff)

    def __str__(self):
        list_monom = list()
        for it in range(0, len(self._coeff), 1):
            tmp_str = GetMonom(it,self._coeff,self.code)
            if tmp_str == "continue":
                continue
            else:
                list_monom.append(tmp_str)
        list_monom.reverse()
        list_monom[0] = EditSeniorMonom(list(list_monom[0]))
        return GetPolinomial(list_monom)

    def __eq__(self, other):
        if type(other) is Polynomial:
            zip(self._coeff, other._coeff)
            for i, j in zip(self._coeff, other._coeff):
                if i != j:
                    return False
            return True

    def __addValue(self, coeff):
        if (type(coeff) is int) or (type(coeff) is float):
            self._coeff.append(coeff)
        elif type(coeff) is Polynomial:
            self._coeff = coeff._coeff
        else:
            self._coeff = list(coeff)
        for it in reversed(self._coeff):
            if it == 0:
                self._coeff.remove(it)
            else:
                break

    def __checkValue(self, coeff):
        #assert (type(coeff) == int)
        if type(coeff) is str:
            raise ValueError("Invalid initialization list! Type Str")
        if type(coeff) is list:
            for it in coeff:
                if (type(it) is not int) and (type(it) is not float):
                    raise ValueError("Error! Invalid type in initialization list")

# ========================================================================================================#
if __name__ == "__main__":
    #str = 'X%' % 2
    #print(str)
    #a = 12654
    #b = map(int, str(a))
    #for i in b:
    #    print(i)
    #p = Polynomial()
    X = Polynomial((0,-1))
    A = Polynomial((1,2.34,-1,0))
    A1 = Polynomial((1.32, 2.34, -5.67, -1, 99, 102, 44, 4.3, 5.5, 66, 101, 103))
    print(X)
    print(A)
    print(A1)
    #print(A1 == A)
    #print('x'+chr(185))
    # На этом упадет!
    #print(repr(Polynom('1')))
    #print(repr(Polynom('1 2 3')))
    #print(repr(Polynom('1 2.34 5.67 0')))
