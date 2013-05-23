function rhs = rightSide(x,y,dummy,m)

% b = 1;
% [ReQ,ImQ] = Qm(x,b);
 rhs = [y(2); -y(2)./x + m*m*y(1)./x./x] ;
% + y(1).*ReQ*m./x./x];
% rhs = [y(2); 0];
