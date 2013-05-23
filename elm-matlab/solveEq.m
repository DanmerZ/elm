m = -11;
xspan = [1 0.1]

A = 15;
init = [1 A];
[x,y] = ode45('rightSide',xspan,init,[],m);

 y(:,1)

plot(x,y(:,1))
