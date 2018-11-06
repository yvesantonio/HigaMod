%-+--------------------------------------------------------------------
% HiModlab is a general purpose Hierarchical model reduction library.
% Copyright (C) 2006-2017  by the HiMod authors (see authors.txt).
%
% This file is part of the HiModlab library.
%
% The HiModlab library is free software: you can use it, redistribute
% it and/or modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation, either
% version 2 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Acknowledgement and Disclaimer
%
% The development of HiModLab was supported by the US National Science
% Foundation (NSF Project, DMS 1419060 ). Any opinions, findings and conclusions
% or recomendations expressed in the source code and documentation are
% those of the authors and do not necessarily reflect the views of the
% National Science Foundation (NSF).
%
%
%
%-+--------------------------------------------------------------------
function [errL2,errH1,aus_figure,solsol,meshx,meshy] = plot_solution_DD_IGA(size_mb,a_ril,b_ril,cutx,hx, u,bc_up,bc_down,Coeff_forma,Dati_geometrici,caso,sol_exx,forma,dforma,fisx,fisy,psJac,psJac2,p,k,soldx,soldy)

persistent calls current_figure
    
if isempty(calls)
    calls = 0;
    current_figure=figure;
else
    figure(current_figure);
end	

if calls < 10
    strcalls = strcat('000',int2str(calls));
elseif calls < 100
    strcalls = strcat('00',int2str(calls));
elseif calls >= 100
    strcalls = strcat('0',int2str(calls));
end

len = sum(psJac);

nd = 1;   %number of domains
colormap jet;

minimo = 0;  massimo = 0;
minY   = 0;  maxY = 0;

sol   = cell(nd,1);
sol_x = cell(nd,1);
sol_y = cell(nd,1);
meshx = cell(nd,1);

for i = 1:nd
 	% ne intervalli e nx estremi
    ne = round((cutx(i+1)-cutx(i))/hx(i));
    nx = ne+1;
    
    nnx = ne*k+p+1-k; % DA AGGIORNARE
    
    % Mesh FEM in x, nodi equispaziati
    meshx{i} = zeros(nnx,1);
    meshx{i}(1) = cutx(i);
    for h=2:nnx
        meshx{i}(h) = meshx{i}(h-1)+hx(i)/2;
    end
    
    M = 40;     % number of poinst for every fiber in y direction
    meshy = linspace(0,1,M);
    
    yvis=length(meshy);
    
	sol{i}   = zeros(2*nnx,yvis);
    sol_x{i} = zeros(2*nnx,yvis);
    sol_y{i} = zeros(2*nnx,yvis);
    
    import Core.BasisHandler
    
    obj_newModalBasis = BasisHandler();
    
    obj_newModalBasis.dimModalBasis = size_mb(i);
    obj_newModalBasis.evalNodesY = meshy;
    obj_newModalBasis.labelUpBoundCond = bc_up{i};
    obj_newModalBasis.labelDownBoundCond = bc_down{i};
    obj_newModalBasis.coeffForm = Coeff_forma;
    
    [mb,mby] = newModalBasis(obj_newModalBasis);
    % The parameter "meshy" is used by the function in order to show
    % in which points the modal base must be evalueted.
    
    [X,Y] = meshgrid(fisx,meshy);

    L_computed = Dati_geometrici.L(X);
    a_computed = Dati_geometrici.a(X);
    
    Y = L_computed.*Y+a_computed;
    
    approx  = [];
    dapprox = [];
    
    % Dovrei farlo per ogni m, verifico se ne vale la pena con il caso 1D
    in        = 0; 
    switch caso
        case 9
            fi  = 1*len;
        case {1,2,3,4,5,6,7,8,10}
            fi  = 1;
    end
    knot      = augknt([in fi],p+1);
    h         = hx;
    nel       = ne;
    ins       = sort(reshape((h:h:1-h)'*ones(1,k),1,k*(nel-1)));
    ins       = (fi-in)*ins+in;
    [cp,knot] = bspkntins(p,in:(fi-in)*1/p:fi,knot,ins);
    xx2 = linspace(0,1*len,2*nnx);
    xx3 = linspace(0,1,2*nnx)';
    
    switch caso
        case 9
            for m = 1:size_mb
                
               approx((m-1)*2*nnx+1:2*m*nnx)  = bspeval(p,u{i}((m-1)*nnx+1:m*nnx)',knot,xx2);
               [dc,dk] = bspderiv(p,u{i}((m-1)*nnx+1:m*nnx)',knot);
               dapprox((m-1)*2*nnx+1:2*m*nnx) = bspeval(p-1,dc,dk,xx2);
               
               for j = 2:2*nnx
                    dapprox((m-1)*2*nnx+j) = dapprox((m-1)*2*nnx+j);
               end
            end;
            u{i}  = approx';
            ux = dapprox';
        case {1,2,3,4,5,6,7,8,10}
            for m = 1:size_mb
               approx((m-1)*2*nnx+1:2*m*nnx)  = bspeval(p,u{i}((m-1)*nnx+1:m*nnx)',knot,xx3);
               [dc,dk] = bspderiv(p,u{i}((m-1)*nnx+1:m*nnx)',knot);
               dapprox((m-1)*2*nnx+1:2*m*nnx) = bspeval(p-1,dc,dk,xx3);
               for j = 2:2*nnx
                    dapprox((m-1)*2*nnx+j) = dapprox((m-1)*2*nnx+j)/sqrt(1+(dforma(xx3(j)))^2);
               end
            end;
            u{i}  = approx';
            ux = dapprox';
    end

    solsol = approx';
    for h=1:2*nnx
        for k=1:M
                for imb=1:size_mb(i)
                    sol{i}(h,k)   = sol{i}(h,k)   + u{i}(h+(imb-1)*2*nnx)*mb(k,imb);
                    sol_x{i}(h,k) = sol_x{i}(h,k) + ux(h+(imb-1)*2*nnx)*mb(k,imb);
                    sol_y{i}(h,k) = sol_y{i}(h,k) + u{i}(h+(imb-1)*2*nnx)*mby(k,imb);
                end
		sol{i}(h,k)=sol{i}(h,k)+ a_ril(i)*meshy(k)+b_ril(i);
        end
    end
    
    save('sol','sol')
    
    max(max(sol{i}))
    min(min(sol{i}))
    
    % Codice che permette di costruire il dominio corretto.
    % Riempimento delle matrici X ed Y entrambe di dimensioni (M,nnx)
    % (Le dimensioni sono opposte perch??? poi la sol ??? trasposta)

    X = zeros(M,2*nnx);
    Y = zeros(M,2*nnx);
    k = linspace(cutx(1),cutx(end),2*nnx);
    L=1;
    inv = 1;

    for kk=1:2*nnx
        % Calcolo elementi riga e colonna
        pp  = k(kk);
        if kk>1
            pp1 = k(kk-1);
            if ((-1/dforma(pp))*(-1/dforma(pp1))<0)
                inv = inv + 1;
            end
        end
        xx = linspace(pp+(-1)^(inv)*0.5*sin(-pi/2+atan(-1/dforma(pp))),pp-(-1)^(inv)*0.5*sin(-pi/2+atan(-1/dforma(pp))),M);
        yy = linspace(forma(pp)-(-1)^(inv)*0.5*cos(-pi/2+atan(-1/dforma(pp))),forma(pp)+(-1)^(inv)*0.5*cos(-pi/2+atan(-1/dforma(pp))),M);
        X(:,kk) = xx;
        Y(:,kk) = yy;
    end

    minX = min(min(min(X),0));
    maxX = max(max(max(X),0));
    minY = min(min(min(Y),0));
    maxY = max(max(max(Y),0));
    %[p,e,t] = initmesh([X,Y]);
    %write_msh(p,e,t,'prova');
    %[X,Y] = meshgrid(linspace(0,1,2*nnx),meshy-0.5);
    figure;
    surf(X,Y,sol{1}')
    axis([minX maxX minY maxY]);
    axis equal
    set(gca, 'FontSize', 14)
    figure;
    minimo  = min(min(min(sol{i})),minimo);
    massimo = max(max(max(sol{i})),massimo);
    %minY = min(min(min(Y)),minY);
    %maxY = max(max(max(Y)),maxY);
    soll = sol{i};
    
    
    %display('Il massimo ??:')
    %disp(massimo)
    %display('Il minimo ??:')
    %disp(minimo)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                            %
    %   Salvataggio dell'immagine per plotMTV    %
    %                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sizeX  = size(X,1);
    sizeY  = size(X,2);
    quad_X = sizeX-1;
    quad_Y = sizeY-1;
    num_Th = 2*(sizeX-1)*(sizeY-1);
    
    info_array_sol = zeros(num_Th*3+num_Th,3);
    count = 1;
    info_sol = sol{1}';
    
    for i_quad = 1 : quad_X
        for j_quad = 1 : quad_Y
            
            % Square 1
            info_array_sol(count,1) = X(i_quad , j_quad);
            info_array_sol(count,2) = Y(i_quad , j_quad);
            info_array_sol(count,3) = -info_sol(i_quad , j_quad);
            count = count + 1;
            
            info_array_sol(count,1) = X(i_quad , j_quad+1);
            info_array_sol(count,2) = Y(i_quad , j_quad+1);
            info_array_sol(count,3) = -info_sol(i_quad , j_quad+1);
            count = count + 1;
            
            info_array_sol(count,1) = X(i_quad+1 , j_quad+1);
            info_array_sol(count,2) = Y(i_quad+1 , j_quad+1);
            info_array_sol(count,3) = -info_sol(i_quad+1 , j_quad+1);
            count = count + 1;
            
            info_array_sol(count,1) = 'a';
            info_array_sol(count,2) = 'a';
            info_array_sol(count,3) = 'a';
            count = count + 1;
            
            % Square 2
            info_array_sol(count,1) = X(i_quad , j_quad);
            info_array_sol(count,2) = Y(i_quad , j_quad);
            info_array_sol(count,3) = -info_sol(i_quad , j_quad);
            count = count + 1;
            
            info_array_sol(count,1) = X(i_quad+1 , j_quad+1);
            info_array_sol(count,2) = Y(i_quad+1 , j_quad+1);
            info_array_sol(count,3) = -info_sol(i_quad+1 , j_quad+1);
            count = count + 1;
            
            info_array_sol(count,1) = X(i_quad+1 , j_quad);
            info_array_sol(count,2) = Y(i_quad+1 , j_quad);
            info_array_sol(count,3) = -info_sol(i_quad+1 , j_quad);
            count = count + 1;
            
            info_array_sol(count,1) = 'a';
            info_array_sol(count,2) = 'a';
            info_array_sol(count,3) = 'a';
            count = count + 1;
        end
    end
    
    save ADR_solution.dat info_array_sol -ascii
    disp('Son passato qui')
end

%WARNING:Cambiati a mano per plot
%minimo = 0.0;
%massimo = 0.16;
%%%%%%%%%%%%%%%%%%%%%
contourLines = linspace(minimo, massimo, 20);
%contourLines = linspace(0.07,0.14, 20);
for i=1:nd
    x=0;
    L_computed = Dati_geometrici.L(x);
    a_computed = Dati_geometrici.a(x);
    meshy=meshy*L_computed+a_computed;
    warning('OFF','MATLAB:contourf:ConstantData');
    %contourf(meshx{i},meshy,sol{i}(:,:)',contourLines);
    warning('OFF','MATLAB:contourf:ConstantData');
	%hold on
end
%figure;
for i=1:nd
    x=0;
    L_computed = Dati_geometrici.L(x);
    a_computed = Dati_geometrici.a(x);
    meshy=meshy*L_computed+a_computed;
    warning('OFF','MATLAB:contourf:ConstantData');
    %subplot(1,3,1);
    %contourf(meshx{i},meshy,sol{i}(:,:)',contourLines);
    warning('OFF','MATLAB:contourf:ConstantData');
	%hold on
end
%hold off
set(gca, ...
    'FontUnits','points', ...
    'FontSize',11);

colorbar;

%axis([cutx(1) cutx(end) minY maxY]);
%axis equal
caxis([minimo massimo]);
%caxis([-5*1e-5,25*1e-5]);
calls = calls + 1;
aus_figure=current_figure;

size_fb  = nnx;
meshx    = linspace(0,cutx(end),size_fb);
meshy    = linspace(-0.5,0.5,M);

%Preparo ad hoc il codice per vedere l'errore.
sol_ex   = zeros(M,2*nnx);
sol_exdx   = zeros(M,2*nnx);
sol_exdy   = zeros(M,2*nnx);
dim      = size_fb;
L=sqrt(2);
L = 0.5*(sqrt(5)+1/2*log(2+sqrt(5)));
meshxmia = [];
meshxmia(1) = 0;


for i = 2 : 2*nnx
    meshxmia(i) = meshxmia(i-1) + psJac2(i-1);%+ sqrt(((fisx(i-1)-fisx(i))^2+(fisy(i-1)-fisy(i))^2));
end
meshymia = linspace(-0.5,0.5,M);

for ii=1:M
    for j=1:2*nnx
        xv          = meshxmia(j); 
        yv          = meshymia(ii);
        sol_ex(ii,j)   = sol_exx(xv,yv);
        sol_exdx(ii,j) = soldx(xv,yv);
        sol_exdy(ii,j) = soldy(xv,yv);
    end
end

minimo  = 0;
massimo = 0;
minimo  = min(min(min(sol_ex)),minimo);
massimo = max(max(max(sol_ex)),massimo);
contourLines = linspace(minimo,massimo,20);


%subplot(1,3,2);
%figure;
%contourf(X,Y,sol_ex);%,contourLines);
%axis equal
%colorbar;


%figure;
sum(sum(sol_ex.^2))
sum(sum(sol{1}.^2))
err1 = abs(sol_ex'   - sol{1}  );
err2 = abs(sol_exdx' - sol_x{1});
err3 = abs(sol_exdy' - sol_y{1});
errL2 = sqrt(sum(sum(err1.^2))/sum(sum(sol_ex.^2)));
errH1 = sqrt(sum(sum(err1.^2))/sum(sum(sol_ex.^2))  +...
             sum(sum(err2.^2))/sum(sum(sol_exdx.^2))+...
             sum(sum(err3.^2))/sum(sum(sol_exdy.^2)));
%figure;

%contourLines = linspace(min(min(err1)),max(max(err1)),20);
%subplot(1,3,3);

%surf(X,Y,err1');%,contourLines);
%axis equal
%colorbar;

% figure
% plot(meshymia, sol_ex(1,:));
% hold on
% plot(meshymia, sol{1}(1,:),'g');
end
