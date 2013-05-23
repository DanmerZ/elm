m = -11;
xspan = [0.1 1];
A = 5.5e-9;
init = [0 A];

[x,y] = ode45('rightSide',xspan,init,[],m);
y(:,1)

plot(x,y(:,1))
