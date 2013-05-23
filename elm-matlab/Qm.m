function [ReQ,ImQ] = Qm(x,b)
%imax = 100000;
%x = 0.0:1.0/imax:1.0;
m = -11; n = 3;
c = 2.997924580e+10;
ap = 61; aR = 168;
B0 = 20000.0;

V0 = 0.0; %4500000.0
% b = 1.0;

q_ = q(x,b);
Flux = flux(x,q_);
DFlux = dflux(x,Flux);
Fm_ = Fm(x,q_,m,n);
ne_ = n_e(Flux);
Te_ = Te(Flux);
Ti_ = Ti(Flux);
E_ = Ef(Flux);
Pe_ = 1.6*ne_.*Te_;
Pi_ = 1.6*ne_.*Ti_;
P_ = Pe_ + Pi_;
dPe = dflux(x,Pe_); dPi = dflux(x,Pi_); dP = dflux(x,P_);
sigma = @(x) (1.2e+17) * (Te_./500.0).^1.5;

Km_ = (dPi.*Ti_./Pi_/(ap/100.0) - E_)/B0/30000.0 + Fm_.*x*ap*V0/(m*aR*c);
A_ = 80*pi*m*m*x.*(dPi+dPe).*(1.0./q_./q_-8.0)/B0/B0;
d1 = c./(4.0*pi.*sigma(Flux).*x*ap);
znam = 1.0./((m.*Km_.*Fm_.^2).^2 + (A_.*d1).^2);

ReQ = znam.*m.*Km_.*Km_.*A_.*Fm_.*Fm_;
ImQ = znam.*Km_.*A_.*A_.*d1;