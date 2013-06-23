#! /usr/bin/env python

import numpy
import pylab as pl
from Profiles import *
from Q import *
from scipy import *
from scipy.integrate import odeint

#plotQm()
#x = numpy.linspace(1.0e-300,1.,5000)
x = arange(1.e-300,1.,1.e-4)
b = 1.0  # q(x,b)
[ReQ,ImQ] = Qm(x,b) 

#y1_0 = 1.0  # y(0) = y1(0)
#y0_0 = 1.2e+2 # y'(0) = y2(0)

#y0 = [y0_0, y1_0] # [derivation(0), function(0)]

#ReQ = ReQ[::-1]
#ImQ = ImQ[::-1]

#def func(y,t):	
#    m = -11  
 #   c1 = -3.0*y[1]/t-(1.0-m*m)*y[0]/t/t      
  #  return [c1,y[0]]  # 1-2nd deriv= , 2-1st deriv=
	                     # d  [ y1 ]   [ y2  ]
						 #--- [ y2 ] = [ -y1 ]
						 # dt [    ]   [     ] 
#t = arange(1.,1.e-300,-1.e-4)

#y = odeint(func,y0,t)

#print ReQ
#numpy.savetxt("data",zip(t, y[:,1]))
#pl.plot(t,y[:,1]) #,t,y[:,1])
#pylab.ylim(-1,1)
#pl.plot(x,ReQ)
#pl.show()
