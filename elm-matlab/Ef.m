function fn = Ef(x)
x1 = x;
d1 = 7.32;
r1 = 0.59;
f1 = 1.042;
be = 7.71;
ce = 37;
a1 = 8;
		
ga1 = -673855.0;
ga2 = 57.3688;
ga3 = 1.4;
k1 = -0.6;

y = 1:length(x);
for i = 1:length(x)
if (x(i) >= 0.94)
   y(i) = 1000*(d1+r1*(x1(i)-f1).*exp(be+ce*r1*(x1(i)-f1)).*exp(-(a1*r1*(x1(i)-f1)).^2));
else
   y(i) = 1000*(ga1*exp(-ga2*(x1(i) - ga3).^2) + k1*x1(i));
end
end
fn = y;