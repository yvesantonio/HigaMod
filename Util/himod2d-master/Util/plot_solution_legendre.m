function [aus_figure,sol,meshx,meshy] = plot_solution_legendre(size_mb,mesh_x, u,Map,bc_up,bc_down)

persistent calls current_figure

if isempty(calls)
    calls = 0;
    current_figure=figure('visible','off');
%     current_figure=figure('visible','on');
else
    current_figure=figure('visible','off');
%       current_figure=figure('visible','on');
end

colormap jet;
cmin = -0.04;
cmax = +0.1;

% ne intervalli e nx estremi
nx = length(mesh_x);
%ne=nx-1;

M = 30;
meshy=linspace(-1,1,M)';
yvis=length(meshy);

sol=zeros(nx,yvis);

mb = new_modal_basis_legendre(size_mb,meshy,bc_up,bc_down);

for h=1:nx
    for k=1:M
        for imb=1:size_mb
            sol(h,k)=sol(h,k)+ u(h+(imb-1)*nx)*mb(k,imb);
        end
    end
end

minimo=min(min(sol));
massimo=max(max(sol));

display(minimo)
display(massimo)

%contourLines = linspace(minimo, massimo, 20);
contourLines = linspace(cmin,cmax, 20);

[X,Yhat]=meshgrid(mesh_x,meshy);
Y=Map.invhaty(X,Yhat);
warning('OFF','MATLAB:contourf:ConstantData');
contourf(X,Y,sol(:,:)',contourLines);
grid on
warning('ON','MATLAB:contourf:ConstantData');
hold off
set(gca, ...
    'FontUnits','points', ...
    'FontSize',11);

axis equal
axis([mesh_x(1) mesh_x(end) -1.5 1.5]);

caxis([cmin,cmax]);
calls = calls + 1;
aus_figure=current_figure;
%export_fig(current_figure, 'filename',['himod_m=',num2str(size_mb)],'-eps','-transparent')
end
