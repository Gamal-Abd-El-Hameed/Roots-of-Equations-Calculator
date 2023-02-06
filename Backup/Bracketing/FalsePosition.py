import math
from abc import ABC, abstractmethod
from sympy import *
from Bracketing.Bracketing import Bracketing
from Gauss.FloatConverter import FloatConverter


class FalsePosition(Bracketing):
    def __init__(self, converter: FloatConverter, function: str, rel_tolerance: float, x_lower: float, x_upper: float,
                 max_iterations: int):
        super().__init__(converter, function, rel_tolerance, x_lower, x_upper, max_iterations)

    def solve(self, i: int) -> float or None:
        x = symbols("x")
        func_lower, func_upper = self._get_f_x()
        if self.converter.convert(func_upper * func_lower) > 0:
            raise ValueError("Method does not converge.")
        self.x_middle = self.converter.convert(self.x_upper - self.converter.convert(
            func_upper * self.converter.convert(self.converter.convert(self.x_upper - self.x_lower) / self.converter.convert(func_upper - func_lower))))
        if(self._check_tolerance()):
            return self.x_middle
        func_middle = self.converter.convert(N(sympify(self.function).subs(x, self.x_middle), maxn=self.converter.precision))
        if i == self.max_iterations - 1:
            return self.x_middle
        if self.converter.convert(func_middle * func_lower) < 0:
            self.x_upper = self.x_middle
            return self.solve(i + 1)
        elif self.converter.convert(func_middle * func_lower) > 0:
            self.x_lower = self.x_middle
            return self.solve(i + 1)
        elif self.converter.convert(func_middle * func_lower) == 0:
            return self.x_middle

