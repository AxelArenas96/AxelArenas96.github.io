figure
subplot(1,3,1)
syms t
x=piecewise(t<0,0,0<=t<3,1,3<t<=4,-1,t>4,0);
fplot(x,[-0.25,4.25],'r')
grid on
title('x(t)')
axis([-1 6 -1 2]);
subplot(1,3,2)
h=piecewise(t<0,0,0<=t<3,1,3<t<=4,-1,t>4,0);
fplot(h,[-0.25,4.25],'r')
grid on
title('x(t)')
axis([-1 6 -1 2]);
t1=-4:0.01:-3;
t2=-3:0.01:-1;
t3=-1:0.01:0;
t4=0:0.01:1;
t5=1:0.01:3;
t6=3:0.01:4;
subplot(1,3,3)
plot(t1,(-t1-4))
hold on
plot(t2,(t2+2))
plot(t3,(3*t3+4))
plot(t4,(-3*t4+4))
plot(t5,(-t5+2))
plot(t6,(t6-4))
grid on
title('Autocorrelación x(t)*x(-t)')
axis([-4.75 4.75 -1.25 4.25]);