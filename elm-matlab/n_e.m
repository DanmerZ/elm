function fn = n_e(x)
deltaX = 0.0435;
x1 = x - deltaX;
nsep = 1.2;
Delta = 0.086;
an0 = 3.4;	 
xmid = 1 - Delta/2.0;		 
c1 = 0;

xped = 1 - Delta;
an1 = 5; %5;
aln1 = 0.8;
aln2 = 2.1;

he = @(x) (x >= 0);

fn = nsep+an0*(tanh(2*(1-xmid)/Delta)-tanh(2.*((x1-c1)-xmid)/Delta)) + ...
he(1 - (x)/xped).*an1.*(1 - ((x)/xped).^aln1).^aln2; 			


