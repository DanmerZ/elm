#! /usr/bin/env python

import numpy
import pylab as pl
from Profiles import *
from Q import *

#plotQm()
x = numpy.linspace(1.0e-300,1.05,5000)
b = 1.0  # q(x,b)
[ReQ,ImQ] = Qm(x,b) 
