import pylab as pl
import numpy as np
import csv
from math import *

def he(x):
	"""Heaviside function"""
	return x >= 0.0

def ne(x):
	"""Electron density, 10^{19}m^{-3}"""
	deltaX = 0.00
	x1 = x - deltaX
	nsep = 1.3
	Delta = 0.082
	an0 = 0.8	 
	xmid = 1.0 - Delta/2.0		 
	c1 = 0.0
	xped = 1.0 - Delta
	an1 = 1.0 
	aln1 = 0.8
	aln2 = 2.1
	
	return nsep + an0*(np.tanh(2.0*(1.0-xmid)/Delta) - \
		np.tanh(2.0*((x1-c1)-xmid)/Delta)) + \
		he(1.0-x/xped)*an1*(abs(1.0-(x/xped))**aln1)**aln2

def Te(x):
	"""Electron temperature, eV"""
	deltaX = 0.00
	x1 = x - deltaX - 0.000
	tsep = 120.0
	Delta = 0.05
	at0 = 560.0			
	xmid = 1.0 - Delta/2.0
	c1 = 0.00 		 
	xped = 1.0 - Delta
	at1 = 4900.
	alt1 = 0.90
	alt2 = 1.50
	
	return tsep + at0*(np.tanh(2.0*(1.0-xmid)/Delta) - \
		np.tanh(2.0*((x1-c1)-xmid)/Delta)) + \
		he(1.0-x/xped)*at1*(abs(1.0-(x/xped))**alt1)**alt2
	
def Ti(x):
	"""Ion temperature, eV"""
	deltaX = 0.00 - 0.036
	x1 = x + deltaX
	tsep = 260.0
	Delta = 0.125
	at0 = 1080.0			
	xmid = 1.0 - Delta/2.0
	c1 = 0.00		 
	xped = 1.0 - Delta
	at1 = 5000.0
	alt1 = 1.5
	alt2 = 1.
	
	return tsep + at0*(np.tanh(2.0*(1.0-xmid)/Delta) - \
		np.tanh(2.0*((x1-c1)-xmid)/Delta)) + \
		he(1.0-x/xped)*at1*(abs(1.0-(x/xped))**alt1)**alt2

def Ef(x):
	"""Electric field, V/m"""
	x1 = x
	d1 = 20.0
	r1 = 1.3
	r2 = 1.3
	f1 = 1.015
	be = 7.78
	ce = 22.0
	a1 = 0.5
	#		
	ga1 = -673855.0
	ga2 = 57.3688
	ga3 = 1.4
	k1 = -0.6
	result = np.zeros(len(x))
	i = 0
	for i in range(0,len(x)):
		#if x1[i] >= 0.94:
		result[i] = 1000.0*(d1+r2*(x1[i]-f1)*exp(be+ce*r1*(x1[i]-f1))*exp(-(a1*r1*(x1[i]-f1))**2))
		#else:
		#result[i] = 1000.0*(ga1*exp(-ga2*(x1[i]-ga3)**2)+k1*x1[i])		
	return result


def readcsv(namef,fun,k,outp,yaxis,start):	
	f = np.genfromtxt(namef,delimiter=',')
	x = np.empty(f.shape[0])
	F = np.empty(f.shape[0])
	for i in range(f.shape[0]):
		if not np.isnan(f[i,0]):
		    x[i] = f[i,0]
		    F[i] = f[i,1]
		    #print x[i],F[i]

	pl.clf()
	pl.plot(x,F,'k.')
	x1 = np.arange(0.0,1.0,1.e-4)
	pl.plot(x1,k*fun(x1))
	pl.xlabel(r'$a_{N}$')
	pl.ylabel(yaxis)
	pl.xlim(start,1.)
	pl.grid(True)
	pl.savefig(outp)
	#pl.show()


readcsv('data/Te60.csv',Te,1.e-3,"graphics/Te60.png",'Te60,keV',0.)
readcsv('data/Ti60.csv',Ti,1.e-3,"graphics/Ti60.png",'Ti60,keV',0.)
readcsv('data/ne_b.csv',ne,1.e-1,"graphics/ne60.png",'ne60,10^20 m^-3',0.)
readcsv('data/E60.csv',Ef,1.e-3,"graphics/E60.png",'E60,kV/m',0.85)


x = np.arange(0.,1.05,1.e-3)

dx = x[1] - x[0]
ne_ = ne(x)
Te_ = Te(x)
Ti_ = Ti(x)
Pe_ = 1.6*ne_*Te_
Pi_ = 1.6*ne_*Ti_
DPe_ = np.gradient(Pe_,dx)
DPi_ = np.gradient(Pi_,dx)
E_ = Ef(x)
ap = 61.
B0 = 20000.

K1_ = (DPi_*Ti_/Pi_/(0.01*ap)-E_) #/B0/30000.0

#pl.plot(x,K1_)
pl.clf()
pl.plot(x,DPe_+DPi_)
pl.xlabel(r'$a_{N}$')
pl.ylabel('Total pressure gradient,Pa')
pl.xlim(0.85,1.)
pl.grid(True)
pl.savefig("graphics/DP60.png")

pl.clf()
pl.plot(x,DPi_*Ti_/Pi_/(0.01*ap),color='black',label='Ti*dpi/dan*...')
pl.plot(x,E_,color='red',label='E')
pl.grid(True)
pl.xlim(0.85,1.)
pl.ylim(-6.e+4,5.e+4)
pl.legend()
pl.xlabel(r'$a_{N}$')
pl.ylabel('E and another term (60)')
pl.savefig("graphics/dp-E.png")
#pl.show()

