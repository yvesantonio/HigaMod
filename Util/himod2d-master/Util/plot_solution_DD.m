function [aus_figure,sol,meshx,meshy] = plot_solution_DD(size_mb,a_ril,b_ril,cutx,hx, u,bc_up,bc_down,Coeff_forma,Dati_geometrici)

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

nd = length(size_mb);



colormap jet;

minimo = 0;  massimo = 0;
minY   = 0;  maxY = 0;

sol=cell(nd,1);
meshx=cell(nd,1);

for i = 1:nd
 
 	% ne intervalli e nx estremi
    ne = round((cutx(i+1)-cutx(i))/hx(i));
    nx = ne+1;
    
    % Mesh FEM in x, nodi equispaziati
    meshx{i} = zeros(nx,1);
    meshx{i}(1) = cutx(i);
    for h=2:nx
        meshx{i}(h) = meshx{i}(h-1)+hx(i);
    end
    
    M = 30;
    meshy=linspace(0,1,M)';
    yvis=length(meshy);
    
	sol{i}=zeros(nx,yvis);
    
    mb = new_modal_basis(size_mb(i),meshy,bc_up{i},bc_down{i},Coeff_forma,Dati_geometrici.L(0));
    
    [X,Y] = meshgrid(meshx{i},meshy);

    L_computed = Dati_geometrici.L(X);
    a_computed = Dati_geometrici.a(X);
    
    Y = L_computed.*Y+a_computed;
    
    for h=1:nx
        for k=1:M
                for imb=1:size_mb(i)
                    sol{i}(h,k)=sol{i}(h,k)+ u{i}(h+(imb-1)*nx)*mb(k,imb);
                end
					sol{i}(h,k)=sol{i}(h,k)+ a_ril(i)*meshy(k)+b_ril(i);
        end
    end
    minimo  = min(min(min(sol{i})),minimo);
    massimo = max(max(max(sol{i})),massimo);
    minY = min(min(min(Y)),minY);
    maxY = max(max(max(Y)),maxY);
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
    contourf(meshx{i},meshy,sol{i}(:,:)',contourLines);
    warning('OFF','MATLAB:contourf:ConstantData');
	hold on
end
hold off
set(gca, ...
    'FontUnits','points', ...
    'FontSize',11);

colorbar;

axis([cutx(1) cutx(end) minY maxY]);
%axis equal
caxis([minimo massimo]);
%caxis([-5*1e-5,25*1e-5]);
calls = calls + 1;
aus_figure=current_figure;
end
