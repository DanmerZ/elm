function rhs = rightSide(x,y,dummy,m)
rhs = [y(2); -y(2)/x + m*m*y(1)/x/x];
% rhs = [y(2); 0];