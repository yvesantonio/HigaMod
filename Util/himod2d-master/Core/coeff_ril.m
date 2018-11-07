function [a,b] = coeff_ril(bc_up,bc_down,dato_up, dato_down,Coeff_forma,L)
%
%	Calcola i coefficienti del rilevamento lineare dei dati al bordo.
%
if nargin < 6
    L = 1;
end
muval=Coeff_forma.mu(0,0)/L;
coeffrobin=Coeff_forma.coeffrobin;
switch [bc_up,bc_down]
    case 'dirdir',   % Dirichlet Dirichlet
        a=dato_up-dato_down;
        b=dato_down;
    case 'dirrob',   % Dirichlet Robin
        if(coeffrobin + muval == 0)
            display('In coeff_ril: necessario rilevamento quadratico: Not Yet')
        end
        b = (muval*dato_up + dato_down)/(muval + coeffrobin);
        a = (dato_up - b);
    case 'robrob',   % Robin Robin
        if(2*muval*coeffrobin+coeffrobin^2 == 0)
            display('In coeff_ril: necessario rilevamento quadratico: Not Yet')
        end
        a = ( dato_up - dato_down )/(coeffrobin + 2*muval);
        b = ( dato_down + muval*a )/coeffrobin;
    case 'robdir',   % Robin Dirichlet
		if(coeffrobin + muval == 0)
            display('In coeff_ril: necessario rilevamento quadratico: Not Yet')
        end
        b = dato_down;
        a = ( dato_up - coeffrobin*dato_down )/( coeffrobin + muval );
	case 'neuneu'
		if( dato_up==0 && dato_down==0)
			a=0;
			b=0;
		else
			display('In coeff_ril: Non homogeneus Neumann BC not available.')
		end
    otherwise,
        a=0;
        b=0;
        disp('In coeff_ril: condizioni di bordo non riconosciute')
end

