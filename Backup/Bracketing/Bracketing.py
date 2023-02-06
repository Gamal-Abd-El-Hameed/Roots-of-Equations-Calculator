from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
from sympy import *

from Gauss.FloatConverter import FloatConverter


class Bracketing(ABC):
    def __init__(self, converter: FloatConverter, function: str, rel_tolerance: float, x_lower: float, x_upper: float, max_iterations: int):
        self.__converter = converter
        self.__function = function
        self.__tolerance = rel_tolerance
        self.__x_lower = x_lower
        self.__x_upper = x_upper
        self.__max_iterations = max_iterations
        self.__x_middle = 0
        self.__x_middle_old = 0

    @property
    def converter(self) -> FloatConverter:
        return self.__converter

    @converter.setter
    def converter(self, converter: FloatConverter) -> None:
        self.__converter = converter

    @property
    def function(self) -> str:
        return self.__function

    @function.setter
    def function(self, function: str) -> None:
        self.__function = function

    @property
    def tolerance(self) -> float:
        return self.__tolerance

    @tolerance.setter
    def tolerance(self, tolerance: float) -> None:
        self.__tolerance = tolerance

    @property
    def x_lower(self) -> float:
        return self.__x_lower

    @x_lower.setter
    def x_lower(self, x_lower: float) -> None:
        self.__x_lower = x_lower

    @property
    def x_upper(self) -> float:
        return self.__x_upper

    @x_upper.setter
    def x_upper(self, x_upper: float) -> None:
        self.__x_upper = x_upper

    @property
    def max_iterations(self) -> int:
        return self.__max_iterations

    @max_iterations.setter
    def max_iterations(self, max_iterations: int) -> None:
        self.__max_iterations = max_iterations

    @property
    def x_middle(self) -> float:
        return self.__x_middle

    @x_middle.setter
    def x_middle(self, x_middle: float) -> None:
        self.__x_middle = x_middle

    @property
    def x_middle_old(self) -> float:
        return self.__x_middle_old

    @x_middle_old.setter
    def x_middle_old(self, x_middle_old: float) -> None:
        self.__x_middle_old = x_middle_old

    def get_plot(self) -> plt.figure:
        x = symbols("x")
        interval = (self.x_upper - self.x_lower) / 20
        x_values = []
        fx_values = []
        iterate = [self.x_lower + interval * i for i in range(0,21)]
        for i in iterate:
            x_values.append(i)
            fx_values.append(self.converter.convert(N(sympify(self.function).subs(x, i), maxn=self.converter.precision)))
        fig_size = (10, 5)
        fig = plt.figure(figsize=fig_size)
        plt.plot(x_values, fx_values)
        plt.axhline(linewidth=2, color='k')
        plt.axvline(linewidth=2, color='k')
        plt.grid()
        plt.xlim([self.x_lower,self.x_upper])
        return fig

    def _get_f_x(self) -> tuple[float, float]:
        x = symbols("x")
        func_lower = self.converter.convert(N(sympify(self.function).subs(x, self.x_lower), maxn=self.converter.precision))
        func_upper = self.converter.convert(N(sympify(self.function).subs(x, self.x_upper), maxn=self.converter.precision))
        return func_lower, func_upper

    def _check_tolerance(self) -> bool:
        approx_rel_error = self.converter.convert(
            abs(self.converter.convert(self.x_middle - self.x_middle_old) / self.x_middle))
        self.x_middle_old = self.x_middle
        if approx_rel_error <= self.tolerance:
            return True
        else:
            return False

    @abstractmethod
    def solve(self, i: int) -> float or None:
        """
        Function to calculate an approximate root of a provided function using a
        bracketing method with required relative tolerance and bounds provided.

        :param i: Current iteration of the function (Starts at 0)
        :param function_in_x: Function to calculate root for in terms of x [f(x)]
        :param x_lower: Lower bound on the x-values of the function.
        :param x_upper: Upper bound on the x-values of the function.
        :param rel_tolerance: The maximum ratio between the absolute approximate error and the current iteration result
         tolerable, if the difference calculated in current iteration is less, the root is returned.
        :return: The calculated approximate root for the provided function.
        """
        pass
