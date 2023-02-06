from GJ.FloatChopper import FloatChopper
from GJ.FloatRounder import FloatRounder
from GJ.FloatConverter import FloatConverter
from FixedPoint.Plotting import plotcurve

from sympy import*

class FixedPointIteration:
    def __init__(self, equation : Function, initial : float, et : float, maxit : int, converter:FloatConverter):
        self.equation=equation
        self.initial = initial
        self.et = et
        self.maxit = maxit
        self.converter = converter
    # used for chopping and rounding
    def __convert(self,x):
        return self.converter.convert(x)
    # used to evaluate the function passed at a certain x
    def __evaluate(self,value,function):
        x=value
        return self.converter.convert(eval(function))
    #get g(x)
    def __getg_x(self):
       #getting constant by doing the operation
       # constant = f'(x)- integration of f'(x)
       x = symbols('x')
       f = self.equation
       derivative = f.diff(Symbol('x'))
       integeration = integrate(derivative, x)
       constant = simplify(f - integeration)
       # if there is no constant , add a dummy constant to both sides of equation f(x) = 0 --> f(x)+c = c
       if constant == 0:
            constant = -10
       # f(x)+c=0 ---> f(x) = -c
       f = simplify(f - constant)
       # simplifying the equation to get g(x)
       f = simplify(f / x)
       f = simplify(-1 * constant / f)
       print("g(x) = {}".format(f))
       return f
    #implementing FixedPoint method
    def fixedPointIteration(self):
     try:
        x=symbols('x')
          # plotting curve
        #if function has one literal (x**2-4,x**3,etc...), we solve for x & plot g(x) & x = answer
        if self.equation.count(x) == 1:
            res = solve(self.equation, x)[0]
            fig = plotcurve(str(res),res)
            return (res),fig
        x0=self.initial
        count=0     #used to count number of iterations
        error=[]    # array used to store error values for comparison
        g_x=str(self.__getg_x())    #getting g(x)
        while(count<=self.maxit):
            res=self.__evaluate(x0,g_x) #getting Xnew
            error.append(abs(res - x0))   #getting approximate  error
            # if more than two iterations , we compare between consecutive error values
            if(count>=2):
                #if error is increasing then function diverges and terminate the program
               if(error[count]>error[count-1] and error[count-1]>error[count-2]):
                    raise ValueError("The method Diverges")
            #if error is less than tolerance, then return the root
            if(abs(error[count])<=self.et):
                fig = plotcurve(g_x,res)
                return res, fig
            #else update X for next iteration and incrementing counter
            x0=res
            count+=1
        raise ValueError ("Need More Iterations!!")
     except ZeroDivisionError:
         # if division by zero is raised it's because one of the results was not in g(x)'s domain
         #ask the user to change the initial value
         raise ValueError("Change the initial guess")
     except TypeError:
            raise ValueError('MathError!!')