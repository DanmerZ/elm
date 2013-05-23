#! /usr/bin/env python

import numpy
from math import *

def he(x):
	"""Heaviside function"""
	return x >= 0.0

def q(x,b):
	"""Safety factor"""
	cq1 = 3.6
	cq2 = 5.6
	return b + cq1*x**cq2

def ne(x):
	"""Electron density, 10^{-19}m^{-3}"""
	deltaX = 0.0435
	x1 = x - deltaX
	nsep = 1.2
	Delta = 0.086
	an0 = 3.4	 
	xmid = 1.0 - Delta/2.0		 
	c1 = 0.0
	xped = 1.0 - Delta
	an1 = 5.0 
	aln1 = 0.8
	aln2 = 2.1
	
	return nsep + an0*(numpy.tanh(2.0*(1.0-xmid)/Delta) - \
		numpy.tanh(2.0*((x1-c1)-xmid)/Delta)) + \
		he(1.0-x/xped)*an1*(abs(1.0-(x/xped))**aln1)**aln2

def Te(x):
	"""Electron temperature, eV"""
	deltaX = 0.0435
	x1 = x - deltaX - 0.008
	tsep = 100.0
	Delta = 0.07
	at0 = 370.0			
	xmid = 1.0 - Delta/2.0
	c1 = -0.029 		 
	xped = 1.0 - Delta
	at1 = 3500.0
	alt1 = 1.30
	alt2 = 1.50
	
	return tsep + at0*(numpy.tanh(2.0*(1.0-xmid)/Delta) - \
		numpy.tanh(2.0*((x1-c1)-xmid)/Delta)) + \
		he(1.0-x/xped)*at1*(abs(1.0-(x/xped))**alt1)**alt2
	
def Ti(x):
	"""Ion temperature, eV"""
	deltaX = 0.0435
	x1 = x + 0.053
	tsep = 100.0
	Delta = 0.115
	at0 = 330.0			
	xmid = 1.0 - Delta/2.0
	c1 = 0.07		 
	xped = 1.0 - Delta
	at1 = 5500.0
	alt1 = 0.80
	alt2 = 2.0
	
	return tsep + at0*(numpy.tanh(2.0*(1.0-xmid)/Delta) - \
		numpy.tanh(2.0*((x1-c1)-xmid)/Delta)) + \
		he(1.0-x/xped)*at1*(abs(1.0-(x/xped))**alt1)**alt2

def Ef(x):
	"""Electric field, V/m"""
	x1 = x
	d1 = 7.32
	r1 = 0.59
	f1 = 1.042
	be = 7.71
	ce = 37.0
	a1 = 8.0		
	ga1 = -673855.0
	ga2 = 57.3688
	ga3 = 1.4
	k1 = -0.6
	result = numpy.zeros(len(x))
	i = 0
	for i in range(0,len(x)):
		if x1[i] >= 0.94:
			result[i] = 1000.0*(d1+r1*(x1[i]-f1)*exp(be+ce*r1* \
			(x1[i]-f1))*exp(-(a1*r1*(x1[i]-f1))**2))
		else:
			result[i] = 1000.0*(ga1*exp(-ga2*(x1[i]-ga3)**2)+k1*x1[i])		
	return result
	
def NFlux(x,q):
	"""Normalized poloidal flux. 
	Argument q: array q(x,b)"""
	q1 = q	
	norm = numpy.trapz(x/q1,x)
	dx = x[1] - x[0]
	s = 0.0
	result = numpy.zeros(len(x))
	for i in range(0,len(x)):
		s += dx*x[i]/q1[i]
		result[i] = s
	return result / norm						

def DFlux(x, flux):
	"""Derivation of normalized poloidal flux.
	Argument flux: array NFlux(x,q)"""
#	result = numpy.zeros(len(x))	
#	for i in range(1,len(x)):
#		result[i] = (flux[i] - flux[i-1]) / (x[i] - x[i-1])
#	return result
	dx = x[1] - x[0]
	return numpy.gradient(flux,dx)
