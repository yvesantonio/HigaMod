addpath ../


%% caso Dir-Rob
bc_up='dir';
dato_up = 1;
bc_down='rob';
dato_down = 3;

global coeffrobin mu 

coeffrobin = 2;

mu=@(x,y) 1;

[a,b] = coeff_ril(bc_up,bc_down,dato_up, dato_down);

l=@(x) a*x+b;

fplot(l,[0,1])

R_0=@(x) -mu(0,0)*a + coeffrobin*x;

R_0(l(0))

%% caso Rob-Rob
bc_up='rob';
dato_up = 9;
bc_down='rob';
dato_down = 2;

global coeffrobin mu 

coeffrobin = 2;

mu=@(x,y) 1;

[a,b] = coeff_ril(bc_up,bc_down,dato_up, dato_down);

l=@(x) a*x+b;

fplot(l,[0,1])

R_0=@(x) -mu(0,0)*a + coeffrobin*x;
R_1=@(x) mu(0,0)*a+coeffrobin*x;
R_0(l(0))
R_1(l(1))
