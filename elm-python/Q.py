#! /usr/bin/env python

import numpy
import pylab as pl
from math import *
from Profiles import *
from matplotlib import rc

def plotP(x,dPe,dPi,dP): 
    #rc('text',usetex=True)    
    #rc('font', family='serif') 
    #rc('xtick',labelsize=19)    
    #rc('ytick',labelsize=19) 
    pl.clf()
    pl.figure(1,figsize=(8,5))
    pl.plot(x,dPe/1000,color='black',linewidth=3.5,linestyle='--', \
    label=r'$\frac{dp_{e}}{da_{N}}$')
    pl.plot(x,dPi/1000,color='black',linewidth=3.5,linestyle=':', \
    label=r'$\frac{dp_{i}}{da_{N}}$')
    pl.plot(x,dP/1000,color='black',linewidth=3.5,linestyle='-', \
    label=r'$\frac{dp}{da_{N}}$')
    pl.xlim(0.9,1.0)
    pl.ylim(-163,0.0)
    pl.yticks([0,-40,-80,-120,-160])
    #pl.legend(loc='lower left',prop={'size':24})
    pl.xlabel(r'$a_{N}$',fontsize=22)
    pl.ylabel('Pressure gradients, kPa',fontsize=22)
    pl.grid(True)
    #pl.show()
    pl.savefig("graphics/dP.png") #,bbox_inches=0
    
def plotEf(x,E):
    #rc('text',usetex=True)    
    #rc('font', family='serif') 
    #rc('xtick',labelsize=19)    
    #rc('ytick',labelsize=19)
    pl.clf()
    pl.plot(x,E/1000.0,color='black',linewidth=3.5,linestyle='-')
    ax = pl.gca()
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',0))
    pl.xlim(0.905,1.0)
    pl.ylim(-16.2,20)
    pl.xlabel(r'$a_{N}$',fontsize=22)
   # pl.ylabel(r'\textbf{\textit{$E_{0a}$}}, kV/m',fontsize=22)
    pl.ylabel('E, kV/m',fontsize=22)

    pl.grid(True)
    
    #pl.show()
    pl.savefig("graphics/E.png")

def Qm(x,b):
    """Qm - plasma current response."""
    m = -11 
    n = 3
    c = 2.997924580e+10  #light speed, cm/sec
    ap = 61              #minor radius, cm
    aR = 168             #major radius, cm
    B0 = 20000.0         #toroidal magnetic field = 2 Tesla
    
    V0 = 0.0             #toroidal rotation velocity, cm/sec
    
    # {Var}_ : variables with underscore are defined as numpy arrays
    q_ = q(x,b)
    Flux_ = NFlux(x,q_)
    DFlux_ = DFlux(x,Flux_)
    Fm_ = abs(m) / q_ - n
    ne_ = ne(Flux_)
    Te_ = Te(Flux_)
    Ti_ = Ti(Flux_)
    E_ = Ef(Flux_)
    plotEf(x,Ef(x))
    Pe_ = 1.6*ne_*Te_    #pressure, Pa
    Pi_ = 1.6*ne_*Ti_
    P_ = Pe_ + Pi_
    dx = x[1] - x[0]
    dPe_ = numpy.gradient(Pe_,dx)
    dPi_ = numpy.gradient(Pi_,dx)
    dP_ = numpy.gradient(P_,dx)
    plotP(x,dPe_, dPi_, dP_)
    
    sigma = 1.2e+17*(Te_/500.0)**1.5  #conductivity, Spitzer
       
    Km_ = (dPi_*Ti_/Pi_/(0.01*ap)-E_)/B0/30000.0 + Fm_*x*ap*V0/(m*aR*c)
    A_ = 80*numpy.pi*m*m*x*dP_*(-8.0+(q_)**(-2)) /B0/B0
    d1 = c/(4.0*numpy.pi*sigma*x*ap)
    znam = 1.0/((m*Km_*Fm_**2)**2 + (A_*d1)**2)
    
    ReQ = znam*m*Km_*Km_*A_*Fm_*Fm_
    ImQ = znam*Km_*A_*A_*d1
    result = (ReQ, ImQ)
    
    plotAll(x,Te_,Ti_,ne_)
    pl.clf()
    pl.xlabel(r'$a_{N}$')
    pl.grid(True)
    pl.ylabel('E-other term')
    pl.plot(x,dPi_*Ti_/Pi_/(0.01*ap),'k',x,E_,'r')
    pl.savefig('graphics/dp-Enorm.png')
    #pl.show()
    
    return result

def plotQm():
    """Plots Qm."""
    x = numpy.linspace(1.0e-300,1.0,5000)
    b = 1.0  # q(x,b)
    [ReQ,ImQ] = Qm(x,b)    
    
    pl.plot(x,ReQ,color='blue', linewidth=2.0,linestyle="-")    
    pl.plot(x,ImQ,color='red', linewidth=2.0,linestyle="-")    
    
    pl.xlim(0.94,0.965)
    #pl.ylim()    
    
    
    
    pl.savefig("graphics/Q.png",bbox_inches=0)
    pl.show()


def plotAll(x,Te_,Ti_,ne_):
	pl.clf()
	pl.plot(x,1.e-3*Te_,label='Te')
	pl.plot(x,1.e-3*Ti_,label='Ti')
	pl.grid(True)
	pl.xlabel(r'$a_{N}$')
	pl.ylabel('T,keV')
	pl.savefig('graphics/T.png')
	pl.clf()
	pl.plot(x,ne_)
	pl.grid(True)
	pl.xlabel(r'$a_{N}$')
	pl.ylabel('ne,10^20 m-3')
	pl.savefig('graphics/ne.png')   
    
    
    
