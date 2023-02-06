# Roots-of-Equations-Calculator
The aim of this project is to compare and analyze the behavior of the different numerical methods used for calculating the roots of equations:
1) Bisection
2) False-Position
3) Fixed point
4) Newton-Raphson
5) Secant Method.

* The program takes as an input the equation, the technique to use and its required parameters (e.g. interval for the bisection method).
* It then accepts a free-text input for a non-linear equation:<br>
a. The equations containing different function: {poly, exp, cos, sin}.<br>
b. The variable used is “x”
* It Plots the function with the boundary functions in case of bisection and false position, g(x) with y = x in case of fixed point, f’(x) in the remaining cases.
* The user chooses any of the previously mentioned methods to solve the given equation via a drop-down list and Parameters, if it applicable for the chosen solving method.
* The user enters the precision (Number of significant figures), EPS and the max number of iterations otherwise default values are used; Default #SFs = System Default, Default Epsilon = 0.00001, and Default Max Iterations = 50.
