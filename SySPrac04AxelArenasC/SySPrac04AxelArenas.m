%%
% <html>
% <IMG src="logoupii.png" width="55" height="50"/>
% <h2> INSTITUTO POLITÉCNICO NACIONAL </h2>
% <h2> UNIDAD PROFESIONAL INTERDISCIPLINARIA EN INGENIERÍA Y TECNOLOGÍAS AVANZADAS </h2>
% <h2> SEÑALES Y SISTEMAS </h2>
% <h1> Práctica 4 'Convolución y Correlación de señales en tiempo continuo' </h1>
% <p> Profesor: Dr. Rafael Martinez Martinez </p>
% <p> Integrantes del equipo: </p>
% <ol>
% <li>Arenas Caldera Axel Jacobo</li>
% <li>Islas Martinez Porfirio Ezequiel</li>
% <li>Angeles Cornejo Maria de los Angeles</li>
% <li>Moreno Pilar Yael Maximiliano</li>
% </ol>
% <p> Grupo 2TV1 </p>
% <h2> Objetivos:
% <ol>
% <li>Conocer métodos básicos de integración</li>
% <li>Manipilación dde instrucciones en MATLAB</li>
% <li>Simular convoluciones y correlaciones de señales contínuas</li>
% </ol>
% </html>




%% _*Introducción*_
% *Métodos numéricos de integración*
% El área definida por la función   $$\int_{a}^{b}f(x)dx$$
% Una recta L1[(a,0), (a, f (0) ) ], la recta L2[ (b,0), (b, f(b)  )]
% y la recta L3[ (a,0), (b,0) ]
% el teorema fundamental del cálculo (en una de sus versiones). Si f es integrable**, y tiene antiderivada
%
% que es $$ \frac{dF(x)}{dx}=f(x) $$ , entonces.
%
% $$\int_{a}^{b}f(x)dx=F(b)-F(a)$$
%
% Ahora, encontrar la antiderivada F a veces resulta ser muy complicado o imposible,
% es por eso que existen las técnicas de integración: 
%
% * Cambio de variable.
% * Integración por partes
% * Por sustitución Trigonometrica 
%
% Pero, ¿quién es $F(x)$ (la antiderivada)? La
% respuesta es que no existe como función elemental; por tanto por los
% métodos tradicionales no se puede hallar la integral de esta función. 
% Obteniendo un polinomio lo suficientemente cercano a la función original, resultará que integrar 
% f(x) o P(x) nos daría un valor numérico muy aproximado, obviamente con un pequeño error. 
%
% La técnica que se utiliza para integrar en este caso es por medio de
% aproximación a la función mediante polinomios lo suficientemente
% parecidos (con el menor error posible) e integrarlos. 
% 
%% _*Fórmulas cerradas de Newton-Cotes*_
%
% * Se piensa que la función a integrar es derivable un determinado número de veces, la derivada es continua.
% * Por teoría de Lagrange, la función se puede expresar como un polinomio de grado n, 
% que va a coincidir con la función en determinados puntos más un término de error.
%
% $$ f(x)=P(x)+\frac{f^{(n+1)}\xi (x)}{(n+1)!}(x-x_{0})(x-x_{1})\cdots
% (x-x_{i}) $$
%
% con $x\in [a,b]$ y $x_{i}\in [a,b]$
%
%
% Si se integra, se tiene la fórmula cerrada de (n+1) puntos de
% Newton-Cotes. (Ver gráfica). Donde por regla general:
%
% * $x_{0}=a$
% * $x_{n}=b$
% * $h=\frac{b-a}{n}$ -- (paso)
% * $x_{i}=x_{0}+ih$ con $i=0,1, \cdots ,n$
%
% Dependiendo el grado n del polinomio, se
% obtiene lo que se conoce como Regla del trapecio (n=1) o primer fórmula cerrada de Newton-Cotes, Regla de Simpson
% (n=2) o segunda fórmula cerrada de Newton-Cotes, y Regla de los tres
% octavos de Simpson (n=3).
%
%% *Regla del trapecio Compuesta*
% Es  una generalización de la regla de trapecio para obtener una mejor aproximación  de la integral
% y consiste en subdividir el intervalo $\left[  {a,b} \right]$ en $n$ subintervalos,  
% todos de la misma longitud $h  = \frac{{b - a}}{n}$.  Se aplica la fórmula $  \int_a^b f\left( x \right)dx \approx \frac{h}{2}\left[ {f\left( a  \right) + 2\mathop \sum \limits_{k = 1}^{n - 1} f\left( {{x_k}} \right) +  f\left( b \right)} \right]$.
%
%% _*Regla de Simpson 1/3*_
% Método de integración para calcular integrales definidas donde se conectan grupos sucesivos de tres
% puntos sobre la curva mediante parábolas de segundo grado. A las fórmulas que resultan de calcular la integral
% bajo estos polinomios se les llama Reglas de Simpson. Se utiliza la fórmula
% $$\int_a^b f\left( x \right)dx \approx \frac{h}{3}\left[ {f\left( a  \right) + 4f\left( {{x_m}} \right) + f\left( b \right)} \right]$
% Con     $h  = \frac{{b - a}}{2}$%
%
%%  _*Regla de Simpson*_
%
% $$ \int_{-1}^{1}e^{-x^{2}}dx= \frac{1}{3} \left [ e^{-1}+4e^{0}+e^{-2} 
% \right ] - \frac{1}{90} \left [ -4e^{-\xi^{2}}(-4\xi^{4}+12\xi^{2}-3)
% \right ] $$ 
%
% $$ = 1.5785+ \frac{2}{45}e^{-\xi^{2}}(-4\xi^{4}+12\xi^{2}-3)) $$ 
%
% con $$ -1 <\xi <1 $$
%

%% _*Regla de los tres octavos de Simpson*_
%
% $$ \int_{-1}^{1}e^{-x^{2}}dx=\frac{3}{8}(\frac{2}{3})\left [ 
% e^{-1}+3e^{-\frac{1}{4}}+3e^{-\frac{1}{4}}+e^{-1} \right ] 
% -\frac{3}{80}(\frac{2}{3})^{5}\left [ -4e^{-\xi^{2}}(-4\xi^{4}+12\xi^{2}
% -3)\right ] $$
%
% $$ = 1.5261+\frac{8}{405}e^{-\xi^{2}}(-4\xi^{4}+12\xi^{2}-3))$$ 
%
% con $$ -1 <\xi <1 $$
%
%% _*Fórmulas cerradas de Newton-Cotes compuestas*_
%
% Son una extensión de las fórmulas cerradas de Newton-Cotes, en las que
% una función se divide en múltiples intervalos y no solo uno, y a cada
% intervalo se le aplica la regla del trapecio, de Simpson, de tres octavos
% de Simpson, etcétera.
%
%% _*Regla compuesta del trapecio*_
%
% $$\int_{a}^{b}f(x)dx=\frac{h}{2}\left [
% f(a)+2\sum_{j=1}^{n-1}f(x_{j})+f(b)\right ]
% -\frac{b-a}{12}h^{2}f^{(2)}(\mu)$ con $\mu\in(a,b)$$
%
% * $n$ - número de subintervalos
% * $h=\frac{b-a}{n}$ - paso
% * $x_{j}=a+jh\;\;j=0,1,2,\cdots,n$
%
%% _*Regla compuesta de Simpson*_
%
% $$\int_{a}^{b}f(x)dx=\frac{h}{3}\left [ f(a)+2\sum_{j=1}^{(n/2)-1}f(x_{2j})
% +4\sum_{j=1}^{n/2}f(x_{2j-1})+f(b) \right ]
% -\left(\frac{b-a}{180}\right)h^{4}f^{(4)}(\mu)$ con $\mu\in(a,b)$$
%
% * $n$ par - número de subintervalos
% * $h=\frac{b-a}{n}$ - paso
% * $x_{j}=a+jh\;\;j=0,1,\cdots,n$
%
%% _*Ejemplo*_
%
% Continuando con el ejemplo anterior, tenemos lo siguiente:
%
%% _*Regla compuesta del trapecio*_
%
% * $n=10$
% * $h=\frac{1-(-1)}{10}=0.2$
% * $x_{0}=-1;\;x_{1}=-0.8;\;x_{2}=-0.6;\cdots ;x_{9}=0.8;\;x_{10}=1$$
%
% $$ \int_{-1}^{1}e^{-x^{2}}dx=\frac{0.2}{2}\left[
% e^{-1}+2\sum_{j=1}^{n-1}e^{-x_{j}^{2}}+e^{-1}\right]-\frac{2}{12}(0.2)^{2}
% \left[-2e^{-\mu^{2}}(1-2\mu^{2})\right] $$
%
% $$ =1.4887+0.0133e^{-\mu^{2}}(1-2\mu^{2});\; -1<\mu<1 $$
%
%% _*Regla compuesta de Simpson*_
%
% * $n=10$
% * $h=\frac{1-(-1)}{10}=0.2$
% * $x_{0}=-1;\;x_{1}=-0.8;\;x_{2}=-0.6;\cdots ;x_{9}=0.8;\;x_{10}=1$$
%
% $$ \int_{-1}^{1}e^{-x^{2}}dx=\frac{0.2}{3}\left[
% e^{-1}+2\sum_{j=1}^{4}e^{-x_{2j}^{2}}+4\sum_{j=1}^{5}e^{-(x_{2j}-1)^{2}}+e^{-1}
% \right]-\frac{2}{180}(0.2)^{4}\left[-4e^{-4^{2}}(-4\mu^{4}+12\mu^{2}-3)\right]
% $$
%
% $$ =1.4936+0.000071e^{-\mu^{2}}(-4\mu^{4}+2\mu^{2}-3);\; -1<\mu<1 $$
%
%% _*Cuadratura gaussiana*_
%
% Las fórmulas cerradas de Newton-Cotes son solo fórmulas de cuadratura,
% estas tenían cierto grado de exactitud, dependiendo del grado del
% polinomio de interpolación de Lagrange.
% En este caso se utilizan polinomios de Legendre para aproximar integrales
% de funciones que satisfagan cierta cantidad de error (el error sea el
% mínimo) esto es lo que se conoce como cuadratura gaussiana.
% Los polinomios de Legendre se definen entre -1 y 1 y son 
%
% $$ \{ 1, x, x^{2}-\frac{1}{3}, x^{3}-\frac{3}{5}x, x^{4}-\frac{6}{7}x^{2}+\frac{3}{35}, \cdots \} $$
%
% Estos polinomios son interesantes, ya que en el intervalo -1 a 1 tienen
% cierta simetría y son ortogonales. Esto nos da el siguiente resultado:
%
% Si tenemos:
% 
% * $P_{n}(x)$ - polinomio de Legendre de grado n
% * $x_{1},x_{2},\cdots,x_{n}$ raíces de $P_{n}(x)$
% * $$ c_{i}=\int_{-1}^{1}\prod_{j=1\neq i}^{n}\frac{x-x_{j}}{x_{i}-x_{j}}dx $$
%
% Entonces:
%
% * $\int_{-1}^{1}P(x)dx=\sum_{i=1}^{n}c_{i}P(x_{i})$ y $P(x)$ con grado
% menor que $2n$.
%
% Resultado que se traduce en la aproximación numérica de una integral de -1 a 1 de
% cualquier función:
%
% $$\int_{-1}^{1}f(x)dx\approx\sum_{i=1}^{n}c_{i}f(x_{i})$$
%
% Los coeficientes $c_{i}$ y las raíces $x_{i}$ de los polinomios de
% Legendre, ya vienen especificados en tablas, y solo es necesario
% sustituir los valores en la fórmula dependiendo de la funciòn que se
% quere integrar y del grado del polinomio de Legendre que se quiera
% trabajar.
%
% Si se quiere integrar sobre cualquier intervalo, se necesita previamente
% realizar una equivalencia entre una integral definida en cualquier
% intervalo y una definida entre -1 y 1, y entonces calcular la
% aproximación por cuadratura gaussiana a la nueva integral definida entre
% -1 y 1. Esta equivalencia es la siguiente y se halla mediante cambios de
% variable:
%
% $$\int_{a}^{b}f(x)dx=\int_{-1}^{1}f\left ( \frac{(b-a)t+b+a}{2}\right )
% \left (\frac{b-a}{2}\right ) dt$$
%
% Con el ejemplo que estábamos tratando, la integración de hace de -1 a 1,
% por lo tanto no se necesita ningún ajuste. La integral quedaría así
% (utilizando las raíces y coeficientes de polinomios de Legendre escritos
% en tablas):
%
% Para n=2:
%
% $$\int_{-1}^{1}e^{-x^{2}}dx\approx
% e^{-(0.5773502692)^{2}}+e^{-(-0.5773502692)^{2}}\approx 1.4331 $$
%
% Para n=3:
%
% $$\int_{-1}^{1}e^{-x^{2}}dx\approx
% 0.5555555556e^{-(0.7745966692)^{2}}+0.8888888889+0.5555555556e^{-(-0.7745966692)^{2}}\approx
% 1.4986 $$
%

%% Ejercicio 1
% a)Analítico
%


figure
subplot(1,3,1)
syms t
x=piecewise(t<0,0,0<t<1,-t+1,1<t<2,t-1,t>2,0);
fplot(x,[-1,3],'r')
grid on
title('x(t)')
axis([-1 3.6 -2 2]);
subplot(1,3,2)
h=piecewise(t<0,0,0<t<1,1,t>1,0);
fplot(h,[-1,3],'r')
grid on
title('h(t)')
axis([-1 3.6 -2 2]);
t1=0:0.01:1;
t2=1:0.01:2;
t3=2:0.01:3;
subplot(1,3,3)
plot(t1,(((-t1.^2)/2)+t1));
hold on
plot(t2,(t2.^2)-(3.*t2)+(5/2));
plot(t3,((-t3.^2)/2)+(2.*t3)-(3/2));
grid on
title('x(t)*h(t)')
axis([-1 4 -0.1 0.6]);
%%
% 
%  b)Calculado por MATLAB
%  
% 

x=@(t)(-t+1).*(t>=0 & t<1)+(t-1).*(t>=1 & t<2);
h=@(t)(1).*(t>=0 & t<1);
convconm(x,h);
%% Ejercicio 2
%%
% 
%  a)Analítico
%  
% 

figure
subplot(1,3,1)
syms t
x=piecewise(t<0,0,0<t<1,t,1<t<2,1,t>2,0);
fplot(x,[-1,3],'r')
grid on
title('x(t)')
axis([-1 3.6 -1 1.5]);
subplot(1,3,2)
h=piecewise(t<1,0,1<t<3,1,t>3,0);
fplot(h,[-1,4],'r')
grid on
title('h(t)')
axis([-1 4 -1 1.5]);
t1=1:0.01:2;
t2=2:0.01:3;
t3=3:0.01:4;
t4=4:0.01:5;
subplot(1,3,3)
plot(t1,(((t1.^2)/2)-t1+(1/2)));
hold on
plot(t2,(t2)-(3/2));
plot(t3,((-t3.^2)/2)+(3.*t3)-(3));
plot(t4, 5-(t4));
grid on
title('x(t)*h(t)')
axis([-1 6 -0.1 1.5]);

%%
% 
%  b)Calculado por MATLAB
% 
% 
x=@(t)(1).*(t>=1 & t<3);
h=@(t)(t).*(t>=0 & t<1)+(1).*(t>=1 & t<2);
convconm(x,h);
%% Ejercicio 3

%%
% 
%  a)Analítico
%  
% 


figure
subplot(1,3,1)
syms t
x=piecewise(t<0,0,0<=t<3,1,3<t<=4,-1,t>4,0);
fplot(x,[-0.25,4.25],'r')
grid on
title('x(t)')
axis([-1 6 -1.5 1.25]);
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

%%
% 
%  b)Calculado por MATLAB
% 
% 


g = @(t)((t>=0&t<=3)-2.*(t>=3&t<=4)+(t>=3&t<=4))    
f = @(t) -(-(t<=0&t>=-3)+2.*(t<=-3&t>=-4)-(t<=-3&t>=-4))
convconm(g,f)

%% Ejercicio 4

%%
% 
%  a)Analítico
%  
% 


figure
subplot(1,3,1)
syms t
x=piecewise(t<0,0,0<t<3,1,3<t<4,-1,t>4,0);
fplot(x,[-2,6],'r')
grid on
title('x(t)')
axis([-1.5 5 -1.5 1.5]);
subplot(1,3,2)
h=piecewise(t<0,0,0<t<3,1,3<t<4,-1,t>4,0);
fplot(h,[-2,6],'r')
grid on
title('h(t)')
axis([-1.5 5 -1.5 1.5]);
t1=-4:0.01:-3;
t2=-3:0.01:-2;
t3=-2:0.01:-1;
t4=-1:0.01:0;
t5=-0:0.01:1;
t6=1:0.01:2;
t7=2:0.01:4;
subplot(1,3,3)
hold on
plot(t1,-t1-4);
plot(t2,t2+2);
plot(t3,3.*t3+6);
plot(t4,-t4+2);
plot(t5,-3.*t5+2 );
plot(t6,-t6);
plot(t7,t7-4);
grid on
title('Autocorrelación')
axis([-5 5 -2.25 3.25 ]);
%%
% 
%  b)Calculado por MATLAB
%  
% 

g = @(t)((t>=0&t<=2)-2.*(t>=2&t<2)-(t>=2&t<=4)) 
f = @(t) -(-(t<=0&t>=-3)+2.*(t<=-3&t>=-4)-(t<=-3&t>=-4))
convconm(g,f)