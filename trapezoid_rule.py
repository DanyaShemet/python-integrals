import math
from cmath import asin

import scipy
from scipy import optimize
from sympy import diff, symbols, cos, sin
from scipy.misc import derivative

def trapezoid_rule(func, a, b, nseg):
    """Правило трапеций
       nseg - число отрезков, на которые разбивается [a;b]"""

    h = 1.0 * (b - a) / nseg
    sum = 0.5 * (func(a) + func(b))
    for i in range(1, nseg):
        sum += func(a + i * h)

    return sum * h

print(trapezoid_rule(lambda x: ((1-x**2)**2) * 2 ** asin(x), 0, 0.5, 10))

м
#
# def f1(x): return ((1-x**2)**2) * 2**asin(x)
#
# pohidna = derivative(f1, 1.0, dx=1.0 , n=1)
#
# def f(x): return -2 * x**2 + 4 * x
# max_x = scipy.optimize.fmin(lambda x: f(x), 0)
# print(max_x)

