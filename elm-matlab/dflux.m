function fn = dflux(x,flux)
y = 1:length(x);
y(1) = 0.0;
flux1 = flux;
for i = 2:length(x)
    y(i) =( flux1(i) - flux1(i-1) ) / ( x(i) - x(i-1) );
end
fn = y;