function plotQmForDifferentB
imax = 100000;
x = 0.0:1.0/imax:1.0;
b = 1.0;

for j = 1:5
[ReQ,ImQ] = Qm(x,b);

colors = ['k','y','c','r','g','b','k.'];

%subplot(2,1,1);
figure(1);
reqPlot = plot(x,ReQ,colors(j)); hold on;
set(reqPlot,'LineWidth', 2.4);
grid on
axis([0.94 0.965 -6600 100])
xlabel('a_{N}'); ylabel('ReQ_{N}');
text(0.9598,-4350,'V_{||0} = 0km/sec')
%title('ReQ_{N}(a_{N})')


%subplot(2,1,2);
figure(2);
imqPlot = plot(x,ImQ,colors(j)); hold on;
set(imqPlot,'LineWidth', 2.4);
grid on
axis([0.94 0.965 -12000 5000])
xlabel('a_{N}'); ylabel('ImQ_{N}');
text(0.9595,-5500,'V_{||0} = 0km/sec');
%title('ImQ_{N}(a_{N})')

b = b - 0.06;

end

figure(1);
l1 = legend('b=1','b=0.94','b=0.88','b=0.82','b=0.76','b=0.70');
set(l1,'Location','SouthEast');
figure(2);
l2 = legend('b=1','b=0.94','b=0.88','b=0.82','b=0.76','b=0.70');
set(l2,'Location','SouthEast');