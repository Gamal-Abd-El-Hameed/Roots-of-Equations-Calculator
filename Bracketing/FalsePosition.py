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
        """
        Method to calculate an approximate root of a provided function using false-position
        method with required relative tolerance and bounds provided.

        :param i: Current iteration of the function (Starts at 0)
        :return: The calculated approximate root for the provided function.
        """
        # declare x as sympy symbol
        x = symbols("x")
        # get f(x_lower), f(x_upper)
        func_lower, func_upper = self._get_f_x()
        # check convergence
        if func_lower * func_upper > 0:
            raise ValueError("Even number of roots in interval or no roots in interval.")
        # calculate x_middle by getting chord intersection with x-axis
        self.x_middle = self.converter.convert(self.x_upper - self.converter.convert(
            func_upper * self.converter.convert(self.converter.convert(self.x_upper - self.x_lower) / self.converter.convert(func_upper - func_lower))))
        # check tolerance and if less than or equal return calculated x_middle
        if self._check_tolerance():
            return self.x_middle
        # get f(x_middle)
        func_middle = self.converter.convert(
            N(sympify(self.function).subs(x, self.x_middle), maxn=self.converter.precision))
        # if max iterations reached, return x middle
        if i == self.max_iterations - 1:
            return self.x_middle
        # if product of f(x_middle), f(x_lower) is negative
        if self.converter.convert(func_middle * func_lower) < 0:
            # x_middle is new x_upper
            self.x_upper = self.x_middle
            # recursively call self
            return self.solve(i + 1)
        # if product of f(x_middle), f(x_lower) is positive
        elif self.converter.convert(func_middle * func_lower) > 0:
            # x_middle is new x_lower
            self.x_lower = self.x_middle
            # recursively call self
            return self.solve(i + 1)
        # if product of f(x_middle), f(x_lower) is zero
        elif self.converter.convert(func_middle * func_lower) == 0:
            # x_middle is exact solution, return x_middle
            return self.x_middle

