function fn = flux(x,q)
% q1 = q;
q1 = 1.0 + 3.6*x.^5.6;
norm = trapz(x,x./q1);
y = 1:length(x);
s = 0.0;
dx = x(2) - x(1);
for i = 1:length(x)
    s = s + dx*x(i)/q1(i);
    y(i) = s;
end
fn = y / norm;
