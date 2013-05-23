function fn = Te(x)
deltaX = 0.0435;
x1 = x - deltaX - 0.008;

tsep = 100;
Delta = 0.07;
at0 = 370;			
xmid = 1 - Delta/2;
c1 = -0.029 ;
		 
xped = 1 - Delta;
at1 = 3500; %3500.0;
alt1 = 1.30;
alt2 = 1.50;

he = @(x) (x >= 0);
fn = tsep+at0*(tanh(2*(1-xmid)/Delta)-tanh(2.*((x1-c1)-xmid)/Delta)) + ...
he(1 - (x)/xped).*at1.*(1 - ((x)/xped).^alt1).^alt2; 