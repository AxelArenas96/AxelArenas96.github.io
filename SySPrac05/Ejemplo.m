% d0=1/2;
% dn=@(n) (1/(n*pi))*sin(n*pi/2);
% t0=-pi;
% tf=pi;
% f=@(t) 0.*(t>=-pi & t<-pi/2)+1.*(t>=-pi/2 & t<=pi/2)+0.*(t>pi/2 & t<=pi);
% armo=15;
% a=-5*pi;
% b=5*pi;
% sfc(t0,tf,dn,d0,f,armo,a,b)

T=3;
dn=1/T;
a=-8;
b=8;
armo=15;
sfc2(dn,armo,a,b,T)
