function plot_solution_DD_BFS(mR,aRilR,bRilR,Omega1DR,hx,uR,bc_up,bc_downR,Coeff_forma,...
                              mL,aRilL,bRilL,~,uL,bc_downL)

colormap jet;

% ne intervalli e nx estremi
ne = round(1/hx);
nx = ne+1;

% Mesh FEM in x, nodi equispaziati
meshxL = zeros(nx,1);
meshxL(1) = 0;
for h=2:nx
    meshxL(h) = meshxL(h-1)+hx;
end
meshxR=meshxL+1;

M = 60;
meshy=linspace(-1,1,M)';
yvis=length(meshy);

sol=zeros(2*nx,yvis);

mbL = new_modal_basis(mL,meshy,bc_up,bc_downL,Coeff_forma,1);
mbR = new_modal_basis(mR,meshy/2+0.5,bc_up,bc_downR,Coeff_forma,2);



for h=1:nx
    for k=1:M
        if meshy(k) >= 0
            for imb=1:mL
                sol(h,k)=sol(h,k)+ uL(h+(imb-1)*nx)*mbL(k,imb);
            end
            sol(h,k)=sol(h,k)+ aRilL*meshy(k)+bRilL;
        else
            sol(h,k)=nan;
        end
    end
end
for h=1:nx
    for k=1:M
        for imb=1:mR
            sol(h+nx,k)=sol(h+nx,k)+ uR(h+(imb-1)*nx)*mbR(k,imb);
        end
        sol(h+nx,k)=sol(h+nx,k)+ aRilR*meshy(k)+bRilR;
    end
end

%WARNING:Cambiati a mano per plot
minimo = -0.2;
massimo = 1.2;
%%%%%%%%%%%%%%%%%%%%
contourLines = linspace(minimo, massimo, 21);

contourf([meshxL;meshxR],meshy,sol',contourLines);

set(gca, ...
    'FontUnits','points', ...
    'FontSize',11);

colorbar;

axis([0 2 -1 1]);
caxis([minimo massimo]);