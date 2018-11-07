function lambda=compute_autovalues(stlv,dstlv,m,bc_up,bc_down,Coeff_forma)
%
%   Calcola i primi m autovalori del problema di Stourm-Liouville associate alle condizioni
%	al bordo scelte.
%   In realtà non usa la derivata, era stata inserita per usare il metodo di 
%	Newton: si può togliere.
%

%toll   = 1e-12;
lambda = zeros( m, 1);

switch [bc_up,bc_down]
    
    case 'dirdir',   % Dirichlet Dirichlet
        for i= 1 : m
            lambda(i) = i*pi;
        end
    case {'dirrob','robdir'},   % Dirichlet Robin
        for i = 1 : m
            a = -pi/2 + i*pi + 1e-2;
            b =  -1e-2 + i*pi + pi/2; % attenzione
            lambda( i ) = fzero(stlv, [a, b]);

        end
%         fplot(stlv,[0,10])
%         hold on
%         plot(lambda,zeros(size(lambda)),'*r')
%         ylim([-20,20])
%         grid on
%         pause
    case 'robrob',   % Robin Robin
        i=1;
		if(Coeff_forma.coeffrobin<=pi/2) % where is mu?
		sv=1e-10;
		a=sv;
		b=pi/2-sv;
			found=false;
				while(~found && a<1)
						if stlv(a)*stlv(b)<0
							[lambda(i)]=fzero(stlv,[a,b]);
							found=true;
						end
					a=a*10;
					b=b*10;
				end			
			i=i+1;
		end
        j=1;
        while(i <= m)
            a = -pi/2 + j*pi + 1e-8;
            b =  -1e-8 + j*pi+pi/2;
            if stlv(a)*stlv(b)>0
                ac=(a+b)/2;
                bc=ac;
            
                lambda(i)=fzero(stlv,[a,bc]);               
                i=i+1;				
                lambda(i)=fzero(stlv,[ac,b]);
            else				
                lambda(i)=fzero(stlv,[a,b]);
            end
            i=i+1;
            j=j+1;
        end 
        fplot(stlv,[0,50])
        hold on
        plot(lambda,zeros(size(lambda)),'*r')
        ylim([-20,20])
        grid on
        pause
    otherwise,
        disp('In compute_autovalues: condizioni di bordo non riconosciute')
end
