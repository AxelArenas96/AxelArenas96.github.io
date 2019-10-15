funcion = input('Ingrese la funcion \n f(x)=','s');
a = 14.999;
b = 15.001;
n = 8; % Número de trapecios
h = (b-a)/n; % Número de intervalos
f = 0;

% Algoritmo Trapecio Compuesto
for k = 1:n-1
    x = a+h*k;
    f = f+eval(funcion);
end

% Evaluación de la función
f = 2*f;
x = a; 
f = f+eval(funcion); 
x = b; 
f = f+eval(funcion);
f = (h/2)*(f);

% Muestra de resultados
fprintf('El valor aproximado es: %10.5f\n\n',f)