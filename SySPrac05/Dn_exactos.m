syms t
armo = 15;
sf = 0.504;
t0 = pi;
x = exp(-t/2); x(1) = (exp(-pi/2) + 1)/2;
for n = -armo:armo
    dn = (1/t0)*int(exp(-t./2)*exp(-j.*n.*2.*t),t,0,pi);
    %fprintf('valor %d:\t %f + %f\n',n,real(dn),imag(dn));
    fprintf('%d\t%f\n',n,abs(dn))
end