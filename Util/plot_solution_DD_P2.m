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
function [err,aus_figure,solsol,meshx,meshy] = plot_solution_DD_P2(size_mb,a_ril,b_ril,cutx,hx, u,bc_up,bc_down,Coeff_forma,Dati_geometrici,caso,sol_exx,forma,dforma,fisx,fisy,psJac,psJac2,p,k)

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

nd = length(size_mb);   %number of domains
colormap jet;

minimo = 0;  massimo = 0;
minY   = 0;  maxY = 0;

sol   = cell(nd,1);
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
    
    meshxx = [0 0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.87 0.95 1];
    %meshx{i}=linspace(0,0.5*(sqrt(5)+1/2*log(2+sqrt(5))),nnx);
    %     meshx{i}=zeros(nx,1);
    %     for h=2:nx
    %         meshx{i}(h) = meshx{i}(h-1)+hx(i);
    %     end [PENSO DA ELIMINARE]

    M = 40;     % number of poinst for every fiber in y direction
    %nnx;
    meshy = linspace(0,1,M);
    % Anche se sembrerebbe più corretto usare come estremi [-0.5,0.5]
    % in realtà la function usa sempre sin(i*pi*x), quindi non si interessa
    % della traslazione.
    yvis=length(meshy);
    
	sol{i}=zeros(2*nnx,yvis);
    
    mb = new_modal_basis(size_mb(i),meshy,bc_up{i},bc_down{i},Coeff_forma);
    % The parameter "meshy" is used by the function in order to show
    % in which points the modal base must be evalueted.
    
    [X,Y] = meshgrid(fisx,meshy);

    L_computed = Dati_geometrici.L(X);
    a_computed = Dati_geometrici.a(X);
    
    Y = L_computed.*Y+a_computed;
    approx = [];
    % Dovrei farlo per ogni m, verifico se ne vale la pena con il caso 1D
     in = 0; fi = 1*len;
     knot      = augknt([in fi],p+1);
     h         = hx;
     nel       = ne;
     ins       = sort(reshape((h:h:1-h)'*ones(1,k),1,k*(nel-1)));
     ins       = (fi-in)*ins+in;
     [cp,knot] = bspkntins(p,in:(fi-in)*1/p:fi,knot,ins);
     xx2 = linspace(0,1*len,2*nnx);
     for m = 1:size_mb
        approx((m-1)*2*nnx+1:2*m*nnx)  = bspeval(p,u{i}((m-1)*nnx+1:m*nnx)',knot,xx2);
     end;
    u{i} = approx';
    solsol = approx';
    for h=1:2*nnx
        for k=1:M
                for imb=1:size_mb(i)
                    %disp(h+(imb-1)*2*nnx)
                    sol{i}(h,k)=sol{i}(h,k)+ u{i}(h+(imb-1)*2*nnx)*mb(k,imb);
                end
					sol{i}(h,k)=sol{i}(h,k)+ a_ril(i)*meshy(k)+b_ril(i);
        end
    end
    
    % Codice che permette di costruire il dominio corretto.
    % Riempimento delle matrici X ed Y entrambe di dimensioni (M,nnx)
    % (Le dimensioni sono opposte perchè poi la sol è trasposta)

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

    %[X,Y] = meshgrid(linspace(0,1,2*nnx),meshy-0.5);
    surf(X,Y,sol{1}')
    axis([minX maxX minY maxY]);
    axis equal
    
    minimo  = min(min(min(sol{i})),minimo);
    massimo = max(max(max(sol{i})),massimo);
    %minY = min(min(min(Y)),minY);
    %maxY = max(max(max(Y)),maxY);
    soll = sol{i};
    %size(soll);
    display('Il massimo è:')
    disp(massimo)
    display('Il minimo è:')
    disp(minimo)
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

% switch caso
%     case 1
%         sol_exx  = @(x,y) (-y.^2+y).*(2-x).^3;
%     case 2
%         sol_exx = @ (x,y) x.*(x-4).*(1-y).*y;  
%     case {3,5}
%         sol_exx = @(x,y) -y.^2+y;
%     case 4
%         sol_exx = @(x,y) y;
%     case 6
%         sol_exx = @(x,y) (-x.^2+x).*sin(2*pi.*y);
% end

%sol_exx = @(x,y) y.^2.*(1-y).^2.*(0.75-y).*x.*(2-x);

%Preparo ad hoc il codice per vedere l'errore.
sol_ex   = zeros(M,2*nnx);
dim      = size_fb;
L=sqrt(2);
L = 0.5*(sqrt(5)+1/2*log(2+sqrt(5)));
meshxmia = [];
meshxmia(1) = 0;
% Però ora la soluzione non è più da valutare su nodi equispaziati. Quindi
% azzarderei che meshmiax non sia corretta. Provo a sistemarla:
% Il problema è che ora i vari 'nnx' nodi sono scalati ma nessun problema!

for i = 2 : 2*nnx
    meshxmia(i) = meshxmia(i-1) + psJac2(i-1);%+ sqrt(((fisx(i-1)-fisx(i))^2+(fisy(i-1)-fisy(i))^2));
end
meshymia = linspace(-0.5,0.5,M);
%meshxmia = linspace(0,L,2*nnx); %HO RITOCCATO QUI al posto di meshxmia

for ii=1:M
    for j=1:2*nnx
        xv          = meshxmia(j);%(2)*(i-1)/(size_fb-1); 
        yv          = meshymia(ii);%(j-1)/29 ;
        %sol_ex(i,j) = (xv^2-sqrt(8)*xv)*sin(2*pi*(yv+0.5));
        %sol_ex(i,j) = (0.2*xv^5-sqrt(2)*xv^2)*sin(2*pi*(yv+0.5));
        sol_ex(ii,j) = sol_exx(xv,yv);
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
err1 = abs(sol_ex'-sol{1});
err = sqrt(sum(sum(err1.^2))/sum(sum(sol_ex.^2)));
figure;
contourLines = linspace(min(min(err1)),max(max(err1)),20);
%subplot(1,3,3);

surf(X,Y,err1');%,contourLines);
axis equal
colorbar;

% figure
% plot(meshymia, sol_ex(1,:));
% hold on
% plot(meshymia, sol{1}(1,:),'g');
end
