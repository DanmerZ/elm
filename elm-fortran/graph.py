from pylab import *

#Bootstrap current
clf()
B = loadtxt("data/Bootstrap.dat")

grid(True)
xlim(0.9,1.)
plot(B[:,0],B[:,1],'k',label='Re')
plot(B[:,0],B[:,2],'k--',label='Im')
xlabel(r'$a_{N}$')
ylabel('Calculate: (R/a)**0.5 ...... ')
legend()
savefig("graphics/Bootstrap.png")

#Pressure perturbation
clf()
B = loadtxt("data/znam_.dat")
grid(True)
xlim(0.93,0.97)
plot(B[:,0],B[:,1],'k',label='Re')
plot(B[:,0],B[:,2],'k--',label='Im')
xlabel(r'$a_{N}$')
ylabel('Pressure perturbation')
legend()
savefig("graphics/perturbation.png")

#Q
clf()
B = loadtxt("data/Q.dat")
grid(True)
xlim(0.93,0.97)
plot(B[:,0],B[:,1],'k',label=r'$ReQ_{m}$')
plot(B[:,0],B[:,2],'k--',label=r'$ImQ_{m}$')
xlabel(r'$a_{N}$')
ylabel(r'$Q_{m}$')
legend()
savefig("graphics/Q.png")