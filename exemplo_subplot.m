% exemplo subplot
x = linspace(0,10);
y1 = sin(x);
y2 = sin(5*x);

figure
subplot(2,1,1);
plot(x,y1)

subplot(2,1,2);
plot(x,y2)

