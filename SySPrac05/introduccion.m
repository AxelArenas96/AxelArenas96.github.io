% Cuando se quiere obtener la serie de Fourier de una funci�n se necesitan ciertos coeficientes, dos en el caso de la serie de Fourier trigonom�trica y uno en el caso de la serie de Fourier trigonom�trica compacta y de la serie de Fourier exponencial compleja. En el caso particular de la serie de Fourier trigonom�trica y la serie de Fourier exponencial compleja debemos de aplicar una integral para obtener dichos coeficientes, esta integral puede realizarse de manera anal�tica o num�rica.
%
% Los m�todos de integraci�n num�ricos permiten obtener una aproximaci�n al resultado real a trav�s de calcular �reas con evaluaciones de la funci�n y otras caracter�sticas que var�an seg�n el m�todo que se utilice. 
%
% De manera que al aplicar un m�todo de integraci�n num�rico a los c�lculos de los coeficientes para las series de Fourier se obtiene una aproximaci�n a los coeficientes reales de la serie de Fourier. Exploraremos dos casos para ello. 
%
% * Coeficientes de Fourier a trav�s de la regla del trapecio compuesta
%
% La regla del trapecio calcula el �rea de un trapecio definido por la funci�n en un intervalo. En la regla del trapecio compuesta se extiende esta idea al subdividir el intervalo de integraci�n en varios subintervalos y aplicar en cada uno de ellos la regla del trapecio. 
%
% De manera que la f�rmula es la siguiente:
%
% $$ D_n = \frac{1}{T_0}\int_{<T_0>}f(t)e^{-jnw_0t}dt = \frac{h}{2T_0}[f(a)e^{-jnw_0(a)}+2\sum_{i=1}^{k-1}f(x_i)e^{-jnw_0(x_i)}+f(b)e^{-jnw_0(b)}]$$
%
% Donde : 
%
% $$ T_0 = b - a = $$ un periodo de la funci�n
%
% $$ w_0 = \frac{2\pi}{T_0} $$
%
% $$ k = $$ n�mero de subintervalos
%
% $$ h = \frac{T_0}{k} $$
%
% $$ x_i = a + ih, i = 1,2,...,k-1 $$
%
% * Coeficientes de Fourier a trav�s de la transformada de Fourier discreta
%
% En este caso, se sustituye la integral por la expresi�n de la transformada de Fourier discreta, la cual utiliza muestras de una funci�n en un periodo.
%
% $$ D_n = \frac{1}{T_0}\int_{<T_0>}f(t)e^{-jnw_0t}dt = \lim_{T \to \infty}\frac{1}{N_0T}\sum_{k=0}^{N_0-1}f(kT)e^{-jnw_0kT}T $$
%
% Donde : 
%
% $$ T_0 = $$ un periodo de la funci�n
%
% $$ w_0 = \frac{2\pi}{T_0} $$
%
% $$ T =  $$ el tiempo entre cada muestra
%
% $$ N_0 = \frac{T_0}{T}$$
%
% Reduciendo la ecuaci�n y tomando el valor de T muy peque�o, obtenemos el siguiente resultado:
%
% $$ D_n = \frac{1}{N_0}\sum_{k=0}^{N_0-1}f(kT)e^{-jn \Omega_0} $$
%
% Donde:  
%
% $$ \Omega_0 = w_0T = \frac{2\pi}{N_o} $$
% 
% El obtener los coeficientes de Fourier de esta manera tiene la ventaja que cada $n+N_0$ coeficientes, los coeficientes se repiten. 
% 
% Por �ltimo, se sugiere que los valores de $N_0$ sean una potencia de 2 
%%
% referencia
%
% https://grupocarman.com/blog/efecto-aliasing/
%
% Lathi, B. P. (Bhagwandas Pannalal)
% Linear systems and signals/B. P. Lathi.�2nd