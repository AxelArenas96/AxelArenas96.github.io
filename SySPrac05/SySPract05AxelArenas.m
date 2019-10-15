%% Practica 5 Series de Fourier en tiempo continuo
%
% Materia: _Señales y Sistemas_
%
% Profesor: _Dr. Rafael Martínez Martínez_
%
% Grupo: _2TV1_
% 
% Alumnos:
%
% _1.- Arenas Caldera Axel Jacobo_
%
% _2.- Cornejo Martinez María de los Angeles_
%
% _3.- Islas Martinez Porfirio Ezequiel_
%
% _4.- Moreno Pilar Yael Maximiliano_
%
%% Introducción
% Cuando se quiere obtener la serie de Fourier de una función se necesitan ciertos coeficientes, dos en el caso de la serie de Fourier trigonométrica y uno en el caso de la serie de Fourier trigonométrica compacta y de la serie de Fourier exponencial compleja. En el caso particular de la serie de Fourier trigonométrica y la serie de Fourier exponencial compleja debemos de aplicar una integral para obtener dichos coeficientes, esta integral puede realizarse de manera analítica o numérica.
%
% Los métodos de integración numéricos permiten obtener una aproximación al resultado real a través de calcular áreas con evaluaciones de la función y otras características que varían según el método que se utilice. 
%
% De manera que al aplicar un método de integración numérico a los cálculos de los coeficientes para las series de Fourier se obtiene una aproximación a los coeficientes reales de la serie de Fourier. Exploraremos dos casos para ello. 
%
% * Coeficientes de Fourier a través de la regla del trapecio compuesta
%
% La regla del trapecio calcula el área de un trapecio definido por la función en un intervalo. En la regla del trapecio compuesta se extiende esta idea al subdividir el intervalo de integración en varios subintervalos y aplicar en cada uno de ellos la regla del trapecio. 
%
% De manera que la fórmula es la siguiente:
%
% $$ D_n = \frac{1}{T_0}\int_{<T_0>}f(t)e^{-jnw_0t}dt = \frac{h}{2T_0}[f(a)e^{-jnw_0(a)}+2\sum_{i=1}^{k-1}f(x_i)e^{-jnw_0(x_i)}+f(b)e^{-jnw_0(b)}]$$
%
% Donde : 
%
% $$ T_0 = b - a = $$ un periodo de la función
%
% $$ w_0 = \frac{2\pi}{T_0} $$
%
% $$ k = $$ número de subintervalos
%
% $$ h = \frac{T_0}{k} $$
%
% $$ x_i = a + ih, i = 1,2,...,k-1 $$
%
% * Coeficientes de Fourier a través de la transformada de Fourier discreta
%
% En este caso, se sustituye la integral por la expresión de la transformada de Fourier discreta, la cual utiliza muestras de una función en un periodo.
%
% $$ D_n = \frac{1}{T_0}\int_{<T_0>}f(t)e^{-jnw_0t}dt = \lim_{T \to \infty}\frac{1}{N_0T}\sum_{k=0}^{N_0-1}f(kT)e^{-jnw_0kT}T $$
%
% Donde : 
%
% $$ T_0 = $$ un periodo de la función
%
% $$ w_0 = \frac{2\pi}{T_0} $$
%
% $$ T =  $$ el tiempo entre cada muestra
%
% $$ N_0 = \frac{T_0}{T}$$
%
% Reduciendo la ecuación y tomando el valor de T muy pequeño, obtenemos el siguiente resultado:
%
% $$ D_n = \frac{1}{N_0}\sum_{k=0}^{N_0-1}f(kT)e^{-jn \Omega_0} $$
%
% Donde:  
%
% $$ \Omega_0 = w_0T = \frac{2\pi}{N_o} $$
% 
% El obtener los coeficientes de Fourier de esta manera tiene la ventaja que cada $n+N_0$ coeficientes, los coeficientes se repiten. 
% 
% Por último, se sugiere que los valores de $N_0$ sean una potencia de 2
%% Ejemplo 6.1
%
% Función:
%
% $x(t)=e^{-\frac{t}{2}}$
%
% con:  $T_0=\pi \qquad \omega_0=\frac{2\pi}{T_0}=2\frac{rad}{s}$
%
% Coeficientes:
%
% $a_0=\frac{1}{\pi}\int_{0}^{\pi}e^{-\frac{t}{2}}\;dt=0.504$
%
% $a_n=\frac{2}{\pi}\int_{0}^{\pi}e^{-\frac{t}{2}}cos(2nt)\;dt=0.504\frac{2}{1+16n^{2}}$
%
% $b_n=\frac{2}{\pi}\int_{0}^{\pi}e^{-\frac{t}{2}}sen(2nt)\;dt=0.504\frac{8n}{1+16n^{2}}$
%
% Serie de Fourier Trigonométrica
%
% $x(t)=0.504[1+\sum_{n=1}^{\infty}\frac{2}{1+16n^{2}}(cos(2nt)+4nsen(2nt))]$
%
% Para el espectro trigonométrico necesitamos obtener los coeficientes de
% la serie de Fourier Trigonométrica Compacta como sigue:
% 
% Coeficientes y ángulos
% 
% $C_0=a_0=0.504$
%
% $C_n=\sqrt{a_n^2+b_n^2}=\sqrt{(0.504\frac{2}{1+16n^{2}})^2+(0.504\frac{8n}{1+16n^{2}})^2}=0.504(\frac{2}{\sqrt{1+16n^2}})$
%
% $\theta_n=tan^{-1}(\frac{-b_n}{a_n})=tan^{-1}(-4n)=-tan^{-1}(4n)$
%
% Serie de Fourier Trigonométrica Compacta
%
% $x(t)=0.504+0.504\sum_{n=1}^{\infty}\frac{2}{\sqrt{1+16n^{2}}}cos(2nt-tan^{-1}(4n))]$
%
% Espectro trigonométrico de Fourier con 4 armónicos
%
% <<Ejemplo_61_4a.jpg>>
%
% Espectro trigonométrico de Fourier con 15 armónicos
%
% <<Ejemplo_61_15a.jpg>>
% 
% Luego obtenemos $D_0$ y $D_n$ como sigue:
% 
% $D_0=a_0=0.504$
%
% $D_n=\frac{1}{\pi}\int_{0}^{\pi}e^{-\frac{t}{2}}e^{-2njt}\;dt=-\frac{2}{\pi}(e^{-\frac{\pi}{2}}-1)\frac{1}{1+4nj}=0.5042\frac{1}{1+4nj}$
% 
% Con estos datos procedemos a utilizar la función del Apéndice A con 4 armónicos y el resultado es el siguiente:
%
% <<Ej_61_4a.jpg>>
%
% y ahora con 15 armónicos:
% 
% <<Ej_61_15a.jpg>>
%
%% Ejemplo 6.2
%
% Función:
%
% $x(t)=\{\begin{array}{c}6t\qquad-\frac{1}{2}<t<\frac{1}{2}\\-6t+6\qquad\frac{1}{2}<t<\frac{3}{2}\end{array}$
%
% con:  $T_0=2 \qquad \omega_0=\frac{2\pi}{T_0}=\pi$
%
% Coeficientes:
%
% $a_0=\frac{1}{2}\int_{-\frac{1}{2}}^{\frac{3}{2}}x(t)\;dt=0$
%
% $a_n=\frac{2}{2}\int_{-\frac{1}{2}}^{\frac{3}{2}}x(t)cos(n\pi t)\;dt=0$
%
% Debido a la simetría de la señal $x(t)$ los coeficientes $a_0=a_n=0$
%
% Luego obtenemos $D_0$ y $D_n$ como sigue:
%
% $D_0=a_0=0$
%
% $D_n=\frac{1}{2}\int_{-\frac{1}{2}}^{\frac{3}{2}}x(t)e^{-n\pi
% tj}\;dt=3sen(\frac{n\pi}{2})[\frac{n\pi -2j}{n^2\pi^2}]+[-e^{-1.5n\pi
% j}+e^{-0.5n\pi j}][\frac{3}{2n\pi j}+\frac{3}{n^2\pi^2}]$
%
% Así la serie de Fourier exponencial compleja queda:
%
% $x(t)=\sum_{-\infty}^{\infty}[3sen(\frac{n\pi}{2})[\frac{n\pi
% -2j}{n^2\pi^2}]+[-e^{-1.5n\pi j}+e^{-0.5n\pi j}][\frac{3}{2n\pi
% j}+\frac{3}{n^2\pi^2}]]e^{n\pi jt}$
% 
% Con estos datos procedemos a utilizar la función del Apéndice A con 4
% armónicos y el resultado es el siguiente: 
%
% <<Ej_62_4a.jpg>>
%
% y ahora con 15 armónicos:
% 
% <<Ej_62_15a.jpg>>
%
%% Ejemplo 6.4
%
% Funcion:
%
% $x(t)=\{\begin{array}{c}0\qquad-\pi<t<-\frac{\pi}{2}\\1\qquad-\frac{\pi}{2}<t<\frac{\pi}{2}\\0\qquad\frac{\pi}{2}<t<\pi\end{array}$
%
% con:  $T_0=2\pi \qquad \omega_0=\frac{2\pi}{T_0}=1$
%
% Coeficientes exponenciales
%
% $D_0=\frac{1}{2}$
%
% $D_n=\frac{1}{n\pi}sen(\frac{n\pi}{2})$
%
% Serie de Fourier exponencial compleja
%
% $x(t)=\frac{1}{2}+\sum_{-\infty}^{\infty}\frac{1}{n\pi}sen(\frac{n\pi}{2})e^{njt}$
%
% Luego con estos datos procedemos a utilizar la función del Apéndice A con 4
% armónicos y el resultado es el siguiente: 
%
% <<Ej_64_4a.jpg>>
%
% y ahora con 15 armónicos:
% 
% <<Ej_64_15a.jpg>>
%
%% Ejemplo 6.5
%
% Función:
%
% $x(t)=e^{-\frac{t}{2}}$
%
% con:  $T_0=\pi \qquad \omega_0=\frac{2\pi}{T_0}=2\frac{rad}{s}$
%
% Coeficientes exponenciales
%
% $D_0=a_0=0.504$
%
% $D_n=\frac{1}{\pi}\int_{0}^{\pi}e^{-\frac{t}{2}}e^{-2njt}\;dt=-\frac{2}{\pi}(e^{-\frac{\pi}{2}}-1)\frac{1}{1+4nj}=0.5042\frac{1}{1+4nj}$
%
% Serie de fourier exponencial compleja
%
% $x(t)=0.504\sum_{-\infty}^{\infty}\frac{0.504}{1+4nj}e^{2njt}$
%
% Luego con estos datos procedemos a utilizar la función del Apéndice A con 4
% armónicos y el resultado es el siguiente: 
%
% <<Ej_61_4a.jpg>>
%
% y ahora con 15 armónicos:
% 
% <<Ej_61_15a.jpg>>
%
%% Ejemplo 6.7
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
%%
% Serie para 4 armónicos
%
% <<Ej_67_4a.jpg>>
%
% Serie para 15 armónicos
%
% <<Ej_67_15a.jpg>>
%
%
%% C6.2
%
%   x = @(t) mod(t+pi/2,2*pi) <= pi;
%   t = linspace (-2*pi, 2*pi,1000);
%   x = @(t) mod(t+pi/2,2*pi) <= pi;
%   sumterms = zeros(16, length(t)); sumterms(1,:) = 1/2;
%   for n = 1:size(sumterms,1)-1;
%   sumterms(n+1,:) = (2/(pi*n)*sin(pi*n/2))*cos(n*t);
%   end
%   x_N = cumsum (sumterms); figure(1); clf; ind = 0;
%   for N = [0,1:2:size(sumterms, 1)-1]
%   ind = ind+1; subplot (3,3,ind);
%   plot (t,x_N(N+1)) 
%   plot (t,x(t), 'k--'); axis ([-2*pi 2*pi -0.2 1.2]);
%   xlabel ('t'); 
%   ylabel (['x_{',num2str(N),'} (t)']);
%   end
%
% <<Ej_C6.2.jpg>>
%% Apéndice A
%
% *Código Serie de Fourier Exponencial Compleja*
%
% <include>sfc.m</include>
% 
%% Referencias
%
% https://grupocarman.com/blog/efecto-aliasing/
%
% Lathi, B. P. (Bhagwandas Pannalal)
% Linear systems and signals/B. P. Lathi.—2nd 
%
%