import math

from sympy import *
from numpy import *
import matplotlib.pyplot as plt

from NewtonRaphson.FloatConverter import FloatConverter


class NewtonRaphson:
    def __init__(self, function, initialGuess, iterations, eps, float_converter: FloatConverter):
        self.function = function
        self.initialGuess = initialGuess
        self.iterations = iterations
        self.eps = eps
        self.converter = float_converter

    def getApproximateRelativeError(self, newValue, oldValue): # method to calculate the relative error
        if newValue == 0:
            return math.inf   # Return that relative error = infinity
        ans = self.converter.convert(abs(self.converter.convert(self.converter.convert(newValue - oldValue) / newValue)) * 100)    # calculate the relative error
        return self.converter.convert(ans)

    def getDerivative(self, function):
        # differentiate the function by x
        derivative = function.diff(Symbol('x'))
        # convert the expression to a function
        derivative = lambdify(Symbol('x'), derivative)
        return derivative

    def solve(self):
        figure = plt.figure()
        function = self.function
        # get the derivative of the inserted function
        functionDerivative = self.getDerivative(function)
        # convert the expression to function
        function = lambdify(Symbol('x'), function)
        currentX = self.initialGuess
        # perform the method for specified number of iterations
        for iteration in range(self.iterations):
            # compute x(i+1) = x(i) - func(xi) / func'(xi) and plot the function with the newX tangent
            try:
                # raise error in case of division by zero or negative root
                if function(currentX) in (nan, oo, -oo, 0) or functionDerivative(currentX) in (nan, oo, -oo):
                    raise ValueError()
                newX = self.converter.convert(currentX - self.converter.convert(function(currentX) / functionDerivative(currentX)))
                self.function_plotter(function, functionDerivative, newX)
            except:
                raise ValueError("Division by Zero Error")
            # Check Required Accuracy
            relativeError = self.getApproximateRelativeError(newX, currentX)
            if relativeError < self.eps:    # end the iterations if we reach maximum relative error < eps
                return newX, figure
            # check if there is divergence
            if abs(newX - currentX) >= 10 ** 5:
                raise ValueError("Error! Diverge!!")
            currentX = newX
        raise ValueError("Exceeds the maximum number of iterations without teaching the required accuracy")

    # method to plot the inserted function with its boundary function
    def function_plotter(self, function, derivative, point):
        # the x-values which will be substituted to plot the function
        range = linspace(point-5, point+5)
        # plot the function with red line
        plt.plot(range, function(range), 'r')
        # plot the tangent using the slope, x-shifting and y-shifting
        plt.plot(range, derivative(point)*(range - point) + function(point), 'b--')
        # plot vertical line to show point location
        plt.plot([point, point], [0, function(point)], 'b--')
        # show the point on the function plot
        plt.scatter(point, function(point))
        plt.axvline(x=0, c="black")
        plt.axhline(y=0, c="black")
        plt.grid()
