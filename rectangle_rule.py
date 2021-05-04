import math
from cmath import asin


def _rectangle_rule(func, a, b, nseg, frac):
    """Обобщённое правило прямоугольников."""
    h = 1.0 * (b - a) / nseg
    sum = 0.0
    xstart = a + frac * h # 0 <= frac <= 1 задаёт долю смещения точки,в которой вычисляется функция, от левого края отрезка dx
    for i in range(nseg):
        sum += func(xstart + i * h)

    return sum*h

def left_rectangle_rule(func, a, b, nseg=1):
    """Правило левых прямоугольников"""
    return _rectangle_rule(func, a, b, nseg, 0.0)

def right_rectangle_rule(func, a, b, nseg=1):
    """Правило правых прямоугольников"""
    return _rectangle_rule(func, a, b, nseg, 1.0)

def midpoint_rectangle_rule(func, a, b, nseg=1):
    """Правило прямоугольников со средней точкой"""
    return _rectangle_rule(func, a, b, nseg, 0.5)


print(midpoint_rectangle_rule(lambda x: ((1-x**2)**2) * 2 ** asin(x), 0, 0.5, 10))
print(right_rectangle_rule(lambda x: ((1-x**2)**2) * 2 ** asin(x), 0, 0.5, 10))
print(left_rectangle_rule(lambda x: ((1-x**2)**2) * 2 ** asin(x), 0, 0.5, 10))

