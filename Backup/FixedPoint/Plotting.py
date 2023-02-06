from pyclbr import Function
import matplotlib.pyplot as pl
import numpy as np
# used to plot  g(x) = x
def plotcurve(equation,point):
    linex=np.linspace(int(point-5),int(point+5)) # generating x-axis
    liney=[]
    for i in linex:
        x=i
        liney.append(eval((equation))) # generating y-axis, wherer y = f(x)
    f = pl.figure(dpi=120)
    pl.title("Solution of g(x) = x where g(x) = " + str(equation))
    pl.plot(linex,liney,label='g(x)=x')
    pl.plot(linex,linex)
    return f

