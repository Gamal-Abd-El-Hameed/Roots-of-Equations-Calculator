
from FixedPoint.FloatConverter import FloatConverter
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
            error.append(abs(abs(res - x0) / res) * 100)   #getting approximate relative error
            # if more than one iteration , we compute the approximate relative error between old and new x
            if(count>=1):
                #if error is increasing then function diverges and terminate the program
                if(error[count]>error[count-1]):
                    return "diverge"
            #if error is less than tolerance, then return the root
            if(abs(error[count])<=self.et):
                fig = plotcurve(g_x,res)
                return res, fig
            #else update X for next iteration and incrementing counter
            x0=res
            count+=1

