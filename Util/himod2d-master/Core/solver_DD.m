function [u,a_ril,b_ril,modalbasis]=solver_DD(cutx,size_mb,hx,bc_up,bc_down,dato_up,dato_down,Dati,Coeff_forma,Dati_geometrici,gamma,export,coupling)
%
%  function [u,a_ril,b_ril]=solver_DD(cutx,size_mb,hx,bc_up,bc_down,dato_up,dato_down,
%                         Dati,Coeff_forma,Dati_geometrici,gamma,export,coupling)
%
%       Risolve un problema ellittico definito con una
%       dinamica dominante in direzione X.
%       Si assume che ci sia una condizione di Dirichlet al bordo di INFLOW
%       e una condizione di Neumann omogenea all'OUTFLOW.
%       Il solutore riceve in ingresso tutti i dati del problema e
%       resistituisce la soluzione.
%       La soluzione del problema viene approssimata tramite un algoritmo
%       di tipo Domain Decomposition (Robin-Robin).
%       In ogni dominio il problema viene risolto tramite un metodo di tipo
%       Hi-Mod con l'implementazione delle basi istruite.
%
%
%  Input:
%
%  cutx: 		Vettore contenente gli estremi dei domini,
%               ad esempio [0,2,3,4] rappresenta i domini (0,2),(2,3),(3,4)
%
%  size_mb: 	Vettore contente la dimensione della base modale in
%               ogni dominio, ad esempio [7,5,3]
%
%  hx:      	Vettore contente il passo della griglia elementi finiti,
%               ad esempio [0.1, 0.1,0.1]
%
%  bc_up,bc_down: 	Variabili contenti le etichette con la natura della
%           condizione di bordo, ad esempio bc_up='rob', bc_down='dir'
%			i valori possibili sono 'dir' e 'rob', per le condizioni di
%           neumann invece, se non ancora previsto non e' difficile
%           l'implementazione.
%
% 			'rob':  mu du/dnu + coeffrobin u = dato
%			'dir':  u=dato
%
%  dato_up,dato_down: 	Variabili contenti i valori della condizione di
%                       bordo.
%
%  Dati: 		Struttura dati contenente le @-functions della condizione
%                di Dirichlet all'inflow e della forzante.
%
%  Coeff_forma:         Struttura dati contenente tutti le @-functions e
%                       le costanti relative alla forma bilineare.
%
%  Dati_geometrici:	Struttura dati contenente i dati geometrici del
%           dominio: attenzione per ora funziona solo nel caso
%			L=1, a=0, psi_x=0.
%
%  gamma:		Struttura dati contenente i due valori (R,L) dei
%           coefficienti della condizione di Robin utilizzata
%			nell'algoritmo di Domain-Decomposition
%
%  Output:
%
%  u:			Struttura dati contente la soluzione.
%
%  Per un esempio dell'utilizzo di solver_DD, si veda lo script example.m
%




%****************************************************************%
%                  INIZIALIZZAZIONE                              %
%****************************************************************%
% Calcolo numero dei sottodomini
nd=length(cutx)-1;
% Calcolo massima dimensione della base modale
max_size_mb = max(size_mb);

% Inizializzazione a 0 delle variabili che conterranno i valori di bordo alle interfacce.
boundary_l = zeros(max_size_mb,nd);
boundary_r = zeros(max_size_mb,nd);

% Inizializzazione delle strutture che conterranno le etichette relative alla natura delle condizioni di bordo all'interfaccia
BC_l=cell(nd,1);
BC_r=cell(nd,1);

% Nel caso dei bordi fisici, le etichette sono note e sono Dirichlet all'inflow e Neumann (omogeneo) all'outflow
BC_l{1}  = 'dir';
BC_r{nd} = 'neu';

switch (coupling)
    case 'RR'
        % Variabili ausiliarie dove esporteremo il valore di interfaccia da passare al
        % sottodominio adiacente.
        rob_l = zeros(max_size_mb,nd);
        rob_r = zeros(max_size_mb,nd);
        rob_r_old=rob_r;
        rob_l_old=rob_l;
    case 'ND'
        dir_l = zeros(max_size_mb,nd);
        dir_r=dir_l;
        neu_r=dir_l;
        neu_l=dir_l;
        % Condizione di interfaccia
        rho_old = zeros(max_size_mb, nd);
        %rho_old(1:size_mb(1),1) = bound_l;
        rho = rho_old;
        phi_old = zeros(max_size_mb, nd);
        %phi_old(1:size_mb(1),1) = bound_l;
        phi = phi_old;
        theta=0.8;
    case 'DN'
        dir_l = zeros(max_size_mb,nd);
        dir_r=dir_l;
        neu_r=dir_l;
        neu_l=dir_l;
        % Condizione di interfaccia
        rho_old = zeros(max_size_mb, nd);
        %rho_old(1:size_mb(1),1) = bound_l;
        rho = rho_old;
        phi_old = zeros(max_size_mb, nd);
        %phi_old(1:size_mb(1),1) = bound_l;
        phi = phi_old;
        theta=0.8;
    case ''
        if length(cutx)>2
            error('it is necessary do decide a coupling')
        end
end


toll  = 1e-4;   % Tolleranza per la convergenza del ciclo Domain-Decomposition
kmax  =  35; 	% Numero massimo di iterazioni del ciclo Domain-Decomposition

% Struttura dati che conterra' la soluzione suddivisa nei diversi sottodomini
u = cell(nd,1);

% Inizializzazione dei sistemi algebrici nei sottodomini
A_i = cell(nd,1);
b_i = cell(nd,1);

% Inizializzazione del vettore che conterra' le basi modali valutate nei nodi di quadratura
modalbasis=cell(nd,1);

% Inizializzazione dei vettori che conterranno i valori del rilevamento (ay+b) in ogni sottodominio
a_ril=zeros(nd,1);
b_ril=zeros(nd,1);

%****************************************************************%
%               ASSEMBLAGGIO SISTEMA LINEARE COMPLETO            %
%****************************************************************%

for i = 1:nd
    % Assemblaggio del sistema lineare tramite la function build_system
    [A_i{i},b_i{i},modalbasis{i},a_ril(i),b_ril(i),yq,wyq] = ...
        build_system(size_mb(i),cutx(i),cutx(i+1),hx(i),bc_up{i}, bc_down{i}, ...
        dato_up{i},dato_down{i},(i==1)-(i==nd),Coeff_forma,Dati,Dati_geometrici,gamma,nd,coupling); % (i==1)-(i==nd) = 1 se ?? il primo dominio, 0 se ?? centrale, -1 se ?? l'ultimo.
end

% per salvare il pattern di sparsit??
% fig=figure;
% spy(A_i{1});
% save_img(fig);

%****************************************************************%
%                  CICLO DOMAIN DECOMPOSITION                    %
%****************************************************************%
for kiter=1:kmax
    
    % fprintf('Inizio ciclo della %d?? iterazione del metodo DD.\n', kiter);
    
    for i = 1:nd
        %****************************************************************%
        %          INTERPOLAZIONE DATI DI INTERFACCIA                    %
        %****************************************************************%
        
        % In questa serie di if vengono chiamate le function interp_interface
        % e interp_interflow, queste hanno il compito di prendere i dati provenienti
        % dai sottodimini adiacenti (rilevamenti inclusi) e proiettarle sulle
        % basi utilizzate nel sottodominio i-esimo.
        % Nel caso si utilizzino le stesse basi in tutti i sottodomini
        % anche se con un diverso numero di modi, e' possibile ottimizzare
        % il codice ed eliminare le proiezioni.
        % Sono state introdotte pensando a condizioni di bordo miste.
        if(i~=1)
            switch(coupling)
                case 'RR'
                    BC_l{i}='rob';
                case 'ND'
                    BC_l{i}='dir';
                case 'DN'
                    BC_l{i}='neu';
                case ''
                    %do nothing
                otherwise
                    display('Error!')
                    break
            end
            [boundary_l( 1 : size_mb(i) , i )] = interp_interface(    ...
                boundary_l(1:size_mb(i-1),i), ... % Boundary contiene gi?? il dato corretto ma sulle basi sbagliate.
                a_ril(i-1)*yq+b_ril(i-1)*ones(size(yq)), a_ril(i)*yq+b_ril(i)*ones(size(yq)), ...
                modalbasis{i-1}, modalbasis{i},   ...
                wyq,BC_l{i});
        else
            if isa(Dati.dato_dir,'function_handle')
                xin=cutx(1);
                boundary_l(1:size_mb(1),1)=interp_inflow( ...
                    a_ril(1)*yq+b_ril(1)*ones(size(yq)) , modalbasis{1}, ...
                    boundary_l(1:size_mb(1),1),wyq,	Dati.dato_dir(Dati_geometrici.a(xin)+Dati_geometrici.L(xin)*yq));
            else
                myZero=@(x) 0*x;
                boundary_l(1:size_mb(1),1)=interp_inflow( ...
                    a_ril(1)*yq+b_ril(1)*ones(size(yq)) , modalbasis{1}, ...
                    boundary_l(1:size_mb(1),1),wyq,	myZero(yq) );
                boundary_l(1:size_mb(1),1)=boundary_l(1:size_mb(1),1)+Dati.dato_dir(1:size_mb(1));
            end
        end
        if(i~=nd)
            switch(coupling)
                case 'RR'
                    BC_r{i}='rob';
                case 'ND'
                    BC_r{i}='neu';
                case 'DN'
                    BC_r{i}='dir';
                case ''
                    
                otherwise
                    display('Error!')
                    break
            end
            [boundary_r( 1 : size_mb(i) , i )] = interp_interface( ...
                boundary_r(1:size_mb(i+1),i), ...% Boundary contiene gi?? il dato corretto ma sulle basi sbagliate.
                a_ril(i+1)*yq+b_ril(i+1)*ones(size(yq)), a_ril(i)*yq+b_ril(i)*ones(size(yq)), ...
                modalbasis{i+1}, modalbasis{i}, ...
                wyq,BC_r{i});
        else
            if isfield(Dati,'dato_neu') 
                if ~isempty(Dati.dato_neu)
                   boundary_r( 1 : size_mb(i) , i )=Dati.dato_neu( 1 : size_mb(i) , i );
                end
            end
        end
        
        BC_l{1}='dir';
        BC_r{nd}='neu';
        
        %****************************************************************%
        %        IMPOSIZIONE CONDIZIONI DI BORDO INFLOW / OUTFLOW        %
        %****************************************************************%
        
        % I gradi di liberta' relativi al dato di dirichlet in INFLOW
        % vengono eliminati dal sistema e portati al secondo membro generando
        % un sistema ridotto
        
        [Arid, brid] = impose_boundary( size_mb(i), cutx(i), cutx(i+1), hx(i), ...
            BC_l{i}, boundary_l(1:size_mb(i),i), ...
            BC_r{i}, boundary_r(1:size_mb(i),i), ...
            A_i{i},b_i{i} ...
            );
        %****************************************************************%
        %               SOLUZIONE SISTEMA LINEARE RIDOTTO                %
        %****************************************************************%
        
        ur=Arid\brid;
        
        %****************************************************************%
        %		        AGGIUNTI DOF DI DIRICHLET (INFLOW)               %
        %            CALCOLO QUANTITA' DA ESPORTARE PER DD               %
        %****************************************************************%
        [ui,dir_l(1:size_mb(i), i),dir_r(1:size_mb(i), i),neu_l(1:size_mb(i), i),neu_r(1:size_mb(i), i),rob_l(1:size_mb(i), i),rob_r(1:size_mb(i),i)] = ...
            extract_solution( size_mb(i), cutx(i), cutx(i+1), hx(i), ...
            BC_l{i}, boundary_l(1:size_mb(i),i), ...
            BC_r{i}, boundary_r(1:size_mb(i),i), ...
            ur, A_i{i}, b_i{i},                 ...
            Coeff_forma, gamma ...
            );
        
        %****************************************************************%
        %                      SALVATAGGIO SOLUZIONE                     %
        %****************************************************************%
        u{i} = ui;
    end
    %****************************************************************%
    %                  	 PLOT                                        %
    %****************************************************************%
    if(export)
        [fig,solmtv,xmtv,ymtv]=plot_solution_DD(size_mb,a_ril,b_ril,cutx,hx, u,bc_up,bc_down,Coeff_forma,Dati_geometrici);        
%         export_py(size_mb,a_ril,b_ril,cutx,hx, u,bc_up,bc_down,Coeff_forma,['k=',num2str(kiter)],Dati_geometrici.L(0),Dati_geometrici.a(0));
    else
        %fig=plot_solution_DD(size_mb,a_ril,b_ril,cutx,hx, u,bc_up,bc_down,Coeff_forma,Dati_geometrici);
    end
    %****************************************************************%
    %                  ASSEGNAMENTO NUOVI DATI                       %
    %****************************************************************%
    switch (coupling)
        case 'RR'
            for i = 1 : nd
                if(i~=nd)
                    boundary_l(:,i+1) = rob_r(:,i);%assegno il mio valore destro al mio dominio a destra per cui sar?? il dato sinistro.
                end
                if(i~=1)
                    boundary_r(:,i-1) = rob_l(:,i);%assegno il mio valore sinistro al mio dominio a sinistra per cui sar?? il dato destro.
                end
            end
            
            %****************************************************************%
            %                  TEST DI CONVERGENZA                           %
            %                           &                                    %
            %             	      AGGIORNAMENTO                              %
            %****************************************************************%
            
            % Considero che il metodo sia andato a convergenza quando tutte le quantita' all'interfaccia non
            % crescono piu'!!
            if (norm(rob_r_old-rob_r,'inf')/norm(rob_r,'inf')+norm(rob_l_old-rob_l,'inf')/norm(rob_l,'inf') )<toll || nd==1
                break
            end
            
            % Aggiornamento variabili di interfaccia
            rob_r_old=rob_r;
            rob_l_old=rob_l;
        case 'ND'
            for i = 1:nd-1
                rho(:,i+1) = theta*dir_r(:,i)+(1-theta)*rho_old(:,i+1);
                phi(:,i)   = theta*neu_l(:,i+1) + (1-theta)*phi_old(:,i); %attention phi_old(:,i+1)? I don't think so
                boundary_l(:,i+1) = rho(:,i+1);
                boundary_r(1:size_mb(i),i) = -phi(1:size_mb(i),i);
            end
            
            % controllo convergenza TODO controllarlo!
            if norm(phi-phi_old,'inf')/norm(phi,'inf') + norm(rho-rho_old,'inf')/norm(rho,'inf')<toll
                break
            end
            
            rho_old = rho;
            phi_old = phi;
        case 'DN'
            for i = 1:nd-1
                rho(:,i+1) = theta*neu_r(:,i  ) + (1-theta)*rho_old(:,i+1);
                phi(:,i)   = theta*dir_l(:,i+1) + (1-theta)*phi_old(:,i);
                boundary_l(:,i+1) = -rho(:,i+1);
                boundary_r(1:size_mb(i),i) = phi(1:size_mb(i),i);
            end
            
            % controllo convergenza TODO controllarlo!
            if norm(phi-phi_old,'inf')/norm(phi,'inf') + norm(rho-rho_old,'inf')/norm(rho,'inf')<toll
                break
            end
            
            rho_old = rho;
            phi_old = phi;
        case ''
            break
    end
end % End ciclo DD

if (kmax==kiter)
    fprintf('Warning: not converging, kmax= %d, kiter= %d\n',kmax,kiter )
else
    if(nd>1)
        fprintf('Reached convergence in %d DD-iteration\n', kiter);
    else
        fprintf('Ended\n');
    end
end

if(export)
    %****************************************************************%
    %                    EXPORT_MTV                                  %
    %****************************************************************%
    %  	fprintf('Export-->MTV\n');
    %  	export_mtv(['himod_m=',num2str(size_mb)],solmtv,xmtv,ymtv);
    %****************************************************************%
    %                   SALVATAGGIO IMMAGINI                         %
    %****************************************************************%
    %  	fprintf('Export-->png\n');
    %  	save_img(fig,['himod_m=',num2str(size_mb)]);
    fprintf('Export-->py\n');
    export_py(size_mb,a_ril,b_ril,cutx,hx, u,bc_up,bc_down,Coeff_forma,['himod_m=',num2str(size_mb)],Dati_geometrici.L(0),Dati_geometrici.a(0));
end

end % End function
