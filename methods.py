import math
from cmath import asin


class Quadrature:
    """Базовые определения для квадратурных формул"""
    __sum = 0.0
    __nseg = 1  # число отрезков разбиения
    __ncalls = 0  # считает число вызовов интегрируемой функции

    def __restart(func, x0, x1, nseg0, reset_calls=True):
        """Обнуление всех счётчиков и аккумуляторов.
           Возвращает интеграл методом трапеций на начальном разбиении"""
        if reset_calls:
            Quadrature.__ncalls = 0
        Quadrature.__nseg = nseg0
        # вычисление суммы для метода трапеций с начальным разбиением на nseg0 отрезков
        Quadrature.__sum = 0.5 * (func(x0) + func(x1))
        dx = 1.0 * (x1 - x0) / nseg0
        for i in range(1, nseg0):
            Quadrature.__sum += func(x0 + i * dx)

        Quadrature.__ncalls += 1 + nseg0
        return Quadrature.__sum * dx

    def __double_nseg(func, x0, x1):
        """Вдвое измельчает разбиение.
           Возвращает интеграл методом трапеций на новом разбиении"""
        nseg = Quadrature.__nseg
        dx = (x1 - x0) / nseg
        x = x0 + 0.5 * dx
        i = 0
        AddedSum = 0.0
        for i in range(nseg):
            AddedSum += func(x + i * dx)

        Quadrature.__sum += AddedSum
        Quadrature.__nseg *= 2
        Quadrature.__ncalls += nseg
        return Quadrature.__sum * 0.5 * dx

    def trapezoid(func, x0, x1, rtol=1e-10, nseg0=1):
        """Интегрирование методом трапеций с заданной точностью.
           rtol - относительная точность,
           nseg0 - число отрезков начального разбиения"""
        ans = Quadrature.__restart(func, x0, x1, nseg0)
        old_ans = 0.0
        err_est = max(1, abs(ans))
        while (err_est > abs(rtol * ans)):
            old_ans = ans
            ans = Quadrature.__double_nseg(func, x0, x1)
            err_est = abs(old_ans - ans)

        print("Total function calls: " + str(Quadrature.__ncalls))
        print("Result: " + str(ans))
        return ans

    def simpson(func, x0, x1, rtol=1.0e-10, nseg0=1):
        """Интегрирование методом парабол с заданной точностью.
           rtol - относительная точность,
           nseg0 - число отрезков начального разбиения"""
        old_trapez_sum = Quadrature.__restart(func, x0, x1, nseg0)
        new_trapez_sum = Quadrature.__double_nseg(func, x0, x1)
        ans = (4 * new_trapez_sum - old_trapez_sum) / 3
        old_ans = 0.0
        err_est = max(1, abs(ans))
        while (err_est > abs(rtol * ans)):
            old_ans = ans
            old_trapez_sum = new_trapez_sum
            new_trapez_sum = Quadrature.__double_nseg(func, x0, x1)
            ans = (4 * new_trapez_sum - old_trapez_sum) / 3
            err_est = abs(old_ans - ans)

        print("Total function calls: " + str(Quadrature.__ncalls))
        print("Result: " + str(ans))
        return ans

    def romberg(func, x0, x1, rtol=1e-10, nseg0=1, maxcol=5, reset_calls=True):
        """Интегрирование методом Ромберга
           nseg0 - начальное число отрезков разбиения
           maxcol - максимальный столбец таблицы"""
        # инициализация таблицы
        Itable = [[Quadrature.__restart(func, x0, x1, nseg0, reset_calls)]]
        i = 0
        maxcol = max(0, maxcol)
        ans = Itable[i][i]
        error_est = max(1, abs(ans))
        while (error_est > abs(rtol * ans)):
            old_ans = ans
            i += 1
            d = 4.0
            ans_col = min(i, maxcol)
            Itable.append([Quadrature.__double_nseg(func, x0, x1)] * (ans_col + 1))
            for j in range(0, ans_col):
                diff = Itable[i][j] - Itable[i - 1][j]
                Itable[i][j + 1] = Itable[i][j] + diff / (d - 1.0)
                d *= 4.0

            ans = Itable[i][ans_col]
            if (maxcol <= 1):  # методы трапеций и парабол обрабатываются отдельно
                error_est = abs(ans - Itable[i - 1][-1])
            elif (i > maxcol):
                error_est = abs(ans - Itable[i][min(i - maxcol - 1, maxcol - 1)])
            else:
                error_est = abs(ans - Itable[i - 1][i - 1])

        print("Total function calls: " + str(Quadrature.__ncalls))
        return ans

## Quadrature.trapezoid(lambda x: ((1-x**2)**2) * 2 ** asin(x), 0, 0.5, rtol=1e-9, nseg0=10)
Quadrature.simpson(lambda x: ((1-x**2)**2) * 2 ** asin(x), 0, 0.5, rtol=1e-9, nseg0=10)


