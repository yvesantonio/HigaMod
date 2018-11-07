function [mb,mb_y]=new_modal_basis(m, yq, bc_up, bc_down, Coeff_forma,L)
%	NEW_MODAL_BASIS    function [mb,mb_y]=new_modal_basis(m, yq, bc_up, bc_down, Coeff_forma)
%	Valuta le basi modali e le loro derivate nei nodi yq, risolvendo il problema di SL.
%   Sul dominio di riferimento
if nargin < 6
    L=1;
    warning('length of fibers (y-direction) not specified --> L=1');
end
muval = Coeff_forma.mu(0,0)/L;
coeffRobin = Coeff_forma.coeffrobin;

mb   = zeros( length(yq), m);
mb_y = zeros( length(yq), m);

%Determiniamo gli autovalori
lambda = zeros(m,1);
switch [bc_up,bc_down]
    case 'dirdir',   % Dirichlet Dirichlet
        for i=1:m
            lambda(i) = i*pi;
        end
        B = @(lambda) 0;
        
    case 'dirrob',   % Dirichlet Robin
        
        stlv  = @(lambda)  tan(lambda) + muval*lambda/coeffRobin;
        dstlv = @(lambda)  1 + ( tan(lambda)).^2 + muval/coeffRobin;
        B     = @(lambda) -tan(lambda);
        
        lambda=compute_autovalues(stlv,dstlv,m,bc_up,bc_down,Coeff_forma);
        
    case 'robrob',   % Robin Robin
        
        stlv  = @(lambda) 2*muval*lambda + tan(lambda).*(coeffRobin - muval.^2*lambda.^2/coeffRobin);
        dstlv = @(lambda) 2*muval + ( 1 + ( tan(lambda) )^2 ).*(coeffRobin - muval.^2*lambda.^2/coeffRobin) + tan(lambda)*(-2*muval*lambda/coeffRobin);
        B     = @(lambda) muval*lambda/coeffRobin;
        
        lambda=compute_autovalues(stlv,dstlv,m,bc_up,bc_down,Coeff_forma);
        
        %xd=linspace(0,32,1e6);
        %plot(xd,stlv(xd));
        %hold on
        %plot(lambda,zeros(size(lambda)),'dr')
        %ylim([-100 200])
        %xlim([0 10])
        %grid on
        %stop
        
    case 'robdir',   % Robin Dirichlet
        
        stlv  = @(lambda) muval*lambda/coeffRobin + tan(lambda);
        dstlv = @(lambda) muval/coeffRobin + ( 1 + ( tan(lambda) )^2 );
        
        B     = @(lambda) 0;
        
        lambda=compute_autovalues(stlv,dstlv,m,bc_up,bc_down,Coeff_forma);
        
    case 'neuneu'
        for i=1:m
            lambda(i) = (i-1)*pi;
        end
        B = @(lambda) 1;
    case 'perper'
        
    otherwise,
        disp('In new_modal_basis: condizioni di bordo non riconosciute o non ancora disponibili')
end

switch [bc_up,bc_down]
    case {'dirdir','dirrob','robrob','robdir'}
        for n = 1 : m
            b         = B( lambda(n) );
            L2norm    = sqrt( ( 1 + b^2 )/2 + ( b^2 - 1 )*sin( 2*lambda(n) )/( 4*lambda(n) ) + b*( sin( lambda(n) ) )^2/lambda(n) );
            %normalizzazione per il coefficiente A
            for i = 1 : length(yq)
                mb  ( i, n) =           sin( lambda(n)*yq(i) )/L2norm +           b*cos( lambda(n)*yq(i) )/L2norm;
                mb_y( i, n) = lambda(n)*cos( lambda(n)*yq(i) )/L2norm - b*lambda(n)*sin( lambda(n)*yq(i) )/L2norm;
            end
        end
    case 'neuneu'
        mb  ( :, 1) = ones ( length(yq), 1 );
        mb_y( :, 1) = zeros( length(yq), 1 );
        for n = 2 : m
            for i = 1 : length(yq)
                mb  ( i, n) =             sqrt(2)*cos( lambda(n)*yq(i) );
                mb_y( i, n) = - lambda(n)*sqrt(2)*sin( lambda(n)*yq(i) );
            end
        end
    case 'perper'
        [mb,mb_y]=new_modal_basis_FOURIER(m, yq);
end
mb=mb/sqrt(L);
mb_y=mb_y/sqrt(L)/L;
