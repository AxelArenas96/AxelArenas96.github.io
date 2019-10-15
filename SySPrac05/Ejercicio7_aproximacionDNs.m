T_0 = pi; N_0 = 256; T = T_0/N_0; t = (0:T:T*(N_0-1))'; M = 10;
x = exp(-t/2); x(1) = (exp(-pi/2) + 1)/2;
%--------

figure(1)
D_n = fft (x)/N_0; n = [-N_0/2:N_0/2-1]';

for a = 1:1:5
mag_dft(a)=abs(fftshift(D_n(a)));
end
for a = 1:1:5
ang_dft(a)=angle(fftshift(D_n(a)));
end
mag_dft
ang_dft
clf; subplot (2, 2, 1); stem(n, abs(fftshift (D_n)),'k');
axis ([-M M -.1 .6]); xlabel('n'); ylabel('|D_n|');
subplot (2, 2, 2); stem(n, angle(fftshift(D_n)),'k');
axis([-M M -pi pi]); xlabel ('n'); ylabel('\angle D n [rad]');
%------

T_0 = pi; N_0 = 256; T = T_0/N_0; W_0=2*pi/T_0; t = (0:T:T*(N_0-1))'; M = 10;

x_0 =@(t) exp(-t/2);
x_1 =@(t) exp(-t/2)*exp(-1j*W_0*t);
x_2 =@(t) exp(-t/2)*exp(-2j*W_0*t);
x_3 =@(t) exp(-t/2)*exp(-3j*W_0*t);
x_4 =@(t) exp(-t/2)*exp(-4j*W_0*t);

%%
% Hallando los D0,...D4 (n positiva) por trapecio compuesto
%%
figure(2)
n=[0:1:4];
D_n=[trap_com(x_0,0,pi,pi/2)/T_0,trap_com(x_1,0,pi,pi/2)/T_0,trap_com(x_2,0,pi,pi/2)/T_0,trap_com(x_3,0,pi,pi/2)/T_0,trap_com(x_4,0,pi,pi/2)/T_0];
mag_trap=abs(D_n)
ang_trap=angle(D_n)

clf; subplot (2, 2, 1); stem(n, abs( (D_n)),'k');
axis ([-M M -.1 .6]); xlabel('n'); ylabel('|D_n|');
subplot (2, 2, 2); stem(n, angle((D_n)),'k');
axis([-M M -pi pi]); xlabel ('n'); ylabel('\angle D n [rad]');
hold on
%%
% De acuerdo a Lathi, en el ejemplo 6.5, que se desarrolló en el R10, la
% serie de Fourier exponencial compleja tiene Dn:
% $\left.\frac{-1}{\pi (1/2+j2n)}e^{-(1/2+j2n)t}\right |_{t=0}^{\pi}$
% Para llegar al resultado más preciso, no se utilizó el resultado final de
% Lathi que indica $\frac{0.504}{1+j4n}$. Sino que desarrollando, se
% encontró: $\frac{-1}{\pi (1/2+j2n)}(e^{-\frac{\pi}{2}}-1)$
%%
for n=2:1:5
   mag_ex(n)=abs((-exp(-pi/2)+1)/(pi*(0.5+2j*(n-1)))); 
end
mag_ex(1)=abs((-exp(-pi/2)+1)/(pi*(0.5)));

for n=2:1:5
   ang_ex(n)=angle((-exp(-pi/2)+1)/(pi*(0.5+2j*(n-1)))); 
end
ang_ex(1)=angle((-exp(-pi/2)+1)/(pi*(0.5)));

T = table(mag_dft',ang_dft',mag_trap',ang_trap',mag_ex',ang_ex'); T(1:5,:);
T.Properties.RowNames = {'D0','D1','D2','D3','D4'};
T.Properties.VariableNames{'Var1'} = 'Abs_DFT';
T.Properties.VariableNames{'Var2'} = 'Ang_DFT';
T.Properties.VariableNames{'Var3'} = 'Abs_Trap_Com';
T.Properties.VariableNames{'Var4'} = 'Ang_Trap_Com';
T.Properties.VariableNames{'Var5'} = 'Abs_Exacto';
T.Properties.VariableNames{'Var6'} = 'Ang_Exacto'

% En su forma rectangular

DFT=mag_dft.*exp(1j.*ang_dft);
TRAP=mag_trap.*exp(1j*ang_trap);
EXAC=mag_ex.*exp(1j*ang_ex);

D= table(DFT',TRAP',EXAC');
D.Properties.RowNames = {'D0','D1','D2','D3','D4'};
D.Properties.VariableNames{'Var1'} = 'Dn_por_DFT';
D.Properties.VariableNames{'Var2'} = 'Dn_por_Trap_comp';
D.Properties.VariableNames{'Var3'} = 'Dn_exacto'

% Comparación
DFT_EXAC=abs(EXAC-DFT);
TRAP_EXAC=abs(EXAC-TRAP);

CMP= table(DFT_EXAC',TRAP_EXAC');
CMP.Properties.RowNames = {'D0','D1','D2','D3','D4'};
CMP.Properties.VariableNames{'Var1'} = 'ERROR_DFT';
CMP.Properties.VariableNames{'Var2'} = 'ERROR_TRAP'

%Por lo tanto se puede ver que el algoritmo DFT es más preciso que el del
%trapecio compuesto, puesto que hubo menor error que el valor exacto.