function fn = Ti(x)
deltaX = 0.0435;

x1 = x + 0.053;	 		 

tsep = 0.1;
Delta = 0.115;
at0 = 0.33;			
xmid = 1 - Delta/2;
c1 = 0.07;
		 
xped = 1 - Delta;
at1 = 5500; %5500.0;
alt1 = 0.80;
alt2 = 2.0;

he = @(x) (x >= 0);
fn = 1000*(tsep+at0*(tanh(2*(1-xmid)/Delta)-tanh(2.*((x1-c1)-xmid)/Delta))) + ...
he(1 - (x)/xped).*at1.*(1 - ((x)/xped).^alt1).^alt2; 
