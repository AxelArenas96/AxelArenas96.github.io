clc;clear;
funcion=input('ingrese la funcion \n f(x)=','s');
 b=input('ingrese el limite superior de la integral\n');
 a=input('ingrese el limite inferior de la integral\n');
 n=input('ingrese el numero de intervalos\n');
 h=(b-a)/n;
 f=0;
 for k=1:n-1
 x=a+h*k;
 f=f+eval(funcion);
 end
 f=2*f;
 x=a; f=f+eval(funcion); x=b; f=f+eval(funcion);
 f=(h/2)*(f);
 fprintf('El valor aproximado es: %10.15f\n\n',f)