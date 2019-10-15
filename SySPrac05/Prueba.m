d=0:0.001:2;
f=@(t) exp(-t.^2);

direc=integral(f,-1,1)
trap=trapecio_compuesto(f,-1,1,0)