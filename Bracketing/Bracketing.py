import math
from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
from sympy import *
from Gauss.FloatConverter import FloatConverter

"""
Abstract Class that defines helper functions for children to
find roots using bracketing methods.
All mathematical operations are converted using the specified float converter
"""


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
        """
        Getter method for the converter property.

        :return: Currently set converter.
        """
        return self.__converter

    @converter.setter
    def converter(self, converter: FloatConverter) -> None:
        """
        Setter method for the converter property.

        :param converter: Float converter to set as property of object.
        :return: None
        """
        self.__converter = converter

    @property
    def function(self) -> str:
        """
        Getter method for the mathematical function property.

        :return: Currently set mathematical function.
        """
        return self.__function

    @function.setter
    def function(self, function: str) -> None:
        """
        Setter method for the mathematical function property.

        :param function: mathematical function as string to set as property of object.
        :return: None
        """
        self.__function = function

    @property
    def tolerance(self) -> float:
        """
        Getter method for the relative tolerance property.

        :return: Currently set relative tolerance.
        """
        return self.__tolerance

    @tolerance.setter
    def tolerance(self, tolerance: float) -> None:
        """
        Setter method for the relative tolerance property.

        :param tolerance: relative tolerance as float to set as property of object.
        :return: None
        """
        self.__tolerance = tolerance

    @property
    def x_lower(self) -> float:
        """
        Getter method for the lower bound x property.

        :return: Currently set lower bound x.
        """
        return self.__x_lower

    @x_lower.setter
    def x_lower(self, x_lower: float) -> None:
        """
        Setter method for the lower bound x property.

        :param x_lower: lower bound x as float to set as property of object.
        :return: None
        """
        self.__x_lower = x_lower

    @property
    def x_upper(self) -> float:
        """
        Getter method for the upper bound x property.

        :return: Currently set upper bound x.
        """
        return self.__x_upper

    @x_upper.setter
    def x_upper(self, x_upper: float) -> None:
        """
        Setter method for the upper bound x property.

        :param x_upper: upper bound x as float to set as property of object.
        :return: None
        """
        self.__x_upper = x_upper

    @property
    def max_iterations(self) -> int:
        """
        Getter method for the max iterations property.

        :return: Currently set max iterations.
        """
        return self.__max_iterations

    @max_iterations.setter
    def max_iterations(self, max_iterations: int) -> None:
        """
        Setter method for the max iterations property.

        :param max_iterations: max iterations as int to set as property of object.
        :return: None
        """
        self.__max_iterations = max_iterations

    @property
    def x_middle(self) -> float:
        """
        Getter method for the calculated x middle property.

        :return: Currently set calculated x middle.
        """
        return self.__x_middle

    @x_middle.setter
    def x_middle(self, x_middle: float) -> None:
        """
        Setter method for the middle x property.

        :param x_middle: middle x as float to set as property of object.
        :return: None
        """
        self.__x_middle = x_middle

    @property
    def x_middle_old(self) -> float:
        """
        Getter method for the previous iteration x middle property.

        :return: Currently set previous iteration x middle.
        """
        return self.__x_middle_old

    @x_middle_old.setter
    def x_middle_old(self, x_middle_old: float) -> None:
        """
        Setter method for the previous iteration x middle property.

        :param x_middle_old: previous iteration x middle as float to set as property of object.
        :return: None
        """
        self.__x_middle_old = x_middle_old

    def get_plot(self) -> plt.figure:
        """
        Method to get plot of function provided from provided lower bound of x
        to provided upper bound of x.

        :return: plt.figure containing plotted function.
        """
        # declare x as sympy symbol
        x = symbols("x")
        # divide distance between upper bound and lower bound into interval
        interval = (self.x_upper - self.x_lower) / 100
        # initialize x values list
        x_values = []
        # initialize f(x) values list
        fx_values = []
        # initialize 101 points to calculate f(x) for and plot.
        iterate = [self.x_lower + interval * i for i in range(0,101)]
        # add points to x values list and f(x) values list
        for i in iterate:
            x_values.append(i)
            fx_values.append(self.converter.convert(N(sympify(self.function).subs(x, i), maxn=self.converter.precision)))
        # initialize plot size
        fig_size = (10, 5)
        # create new figure using initialized plot size
        fig = plt.figure(figsize=fig_size)
        # plot curve using points calculated
        plt.plot(x_values, fx_values)
        # plot x axis
        plt.axhline(linewidth=2, color='k')
        # plot y axis
        plt.axvline(linewidth=2, color='k')
        # plot grid
        plt.grid()
        # limit x axis shown to x values between lower and upper bound
        plt.xlim([self.x_lower, self.x_upper])
        # return plot in figure
        return fig

    def _get_f_x(self) -> tuple[float, float]:
        """
        Method to calculate and return f(x_lower) and f(x_upper)

        :return: tuple containing f(x_lower), f(x_upper) in that order.
        """
        # declare x as sympy symbol
        x = symbols("x")
        # get f(x_lower)
        func_lower = self.converter.convert(N(sympify(self.function).subs(x, self.x_lower), maxn=self.converter.precision))
        # get f(x_upper)
        func_upper = self.converter.convert(N(sympify(self.function).subs(x, self.x_upper), maxn=self.converter.precision))
        # return f(x_lower),f(x_upper)
        return func_lower, func_upper

    def _check_tolerance(self) -> bool:
        """
        Method to check if approximate relative error of current iteration
        is less than or equal provided relative tolerance.

        :return: bool value, true if error less than tolerance, false if not less than tolerance.
        """

        if self.x_middle == 0:
            # if x_middle is zero, consider approximate relative error = infinity to avoid exception.
            approx_rel_error = math.inf
        else:
            # get approximate relative error by dividing modulus approximate absolute error by modulus current x middle
            approx_rel_error = self.converter.convert(
            abs(self.converter.convert(self.x_middle - self.x_middle_old) / self.x_middle))
        # set previous iteration x middle to current iteration x middle
        self.x_middle_old = self.x_middle
        # return true if less than or equal tolerance
        if approx_rel_error <= self.tolerance:
            return True
        # return false otherwise
        else:
            return False

    @abstractmethod
    def solve(self, i: int) -> float or None:
        """
        Method to calculate an approximate root of a provided function using a
        bracketing method with required relative tolerance and bounds provided.

        :param i: Current iteration of the function (Starts at 0)
        :return: The calculated approximate root for the provided function.
        """
        pass
