import math

from sympy import *
from numpy import *
import matplotlib.pyplot as plt

from Secant.FloatConverter import FloatConverter


class Secant:
    def __init__(self, function, oldX, x, iterations, eps, float_converter: FloatConverter):
        self.function = function
        self.ITERATIONS = iterations
        self.EPSILON = eps
        self.converter = float_converter
        self.oldX = oldX
        self.x = x
        self.solve()

    def getAbsoluteRelativeError(self, newValue, oldValue):
        if (newValue == oldValue):
            return 0
        if (newValue == 0):
            return math.inf
        ans = self.converter.convert(abs(self.converter.convert(self.converter.convert(
            newValue - oldValue) / newValue)) * 100)
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
        # y[i - 1] = f(x[i - 1])
        self.oldY = function(self.oldX) 
        # y[i] = f(x[i])
        self.y = function(self.x) 

        for iteration in range(self.ITERATIONS):
            try:
                # x[i + 1] = x[i] - f(x[i]) * (x[i] - x[i - 1]) / (f(x[i]) - f(x[i - 1]))
                self.newX = self.converter.convert(self.x - self.converter.convert(self.y * self.converter.convert(
                    self.converter.convert(self.x - self.oldX) / self.converter.convert(self.y - self.oldY))))
                    # plot the function & the Derivative
                self.function_plotter(function, functionDerivative, self.newX)                
            except:
                raise ValueError("Division by Zero Error")

            # Check Required Accuracy
            self.eps = self.getAbsoluteRelativeError(self.newX, self.x)                
            if(self.eps <= self.EPSILON):
                break
            
            # update (shift) the values of x & y                                   
            self.oldX = self.x                      # x[i - 1] = x[i]
            self.oldY = self.y                      # y[i - 1] = y[i]
            self.x = self.newX                      # x[i] = x[i + 1]
            self.y = function(self.newX)       # y[i] = y[i + 1] 

        return self.newX, figure


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

