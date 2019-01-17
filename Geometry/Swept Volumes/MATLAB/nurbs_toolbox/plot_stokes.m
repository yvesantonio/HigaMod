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
function   plot_stokes(sol,x,Y,ModalBasisX,ModalBasisY,ModalBasisP,out,t,etadown,etaup,mx,my,mp)

persistent fig_cont fig_profileVelocity fig_profile fig_pressure fig_domain
persistent videopressure videocont videoprofile videodomain videoprofileVelocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this part we evaluate the solution on a the qudrature grid

nx=length(x);
ne=nx-1;
ndof=ne+nx;
ny=length(Y(:,1));

Ux=zeros(ndof,ny);
Uy=Ux;
P=zeros(nx,ny);

%X- velocity
for k=1:mx
    MB=ModalBasisX(:,k);
    for j=1:ny
        for i=1:ndof
            Ux(i,j) = Ux(i,j) + sol( i + (k-1)*ndof )*MB(j);
        end
    end
end
display(['min Ux =', num2str(min(min(Ux))),' max Ux =', num2str(max(max(Ux)))] )

%Y- velocity
for k=1:my
    MB=ModalBasisY(:,k);
    for j=1:ny
        for i=1:ndof
            Uy(i,j) = Uy(i,j) + sol( ndof*mx + i + (k-1)*ndof )*MB(j);
        end
    end
end
display(['min Uy =', num2str(min(min(Uy))),' max Uy =', num2str(max(max(Uy)))] )

%Pressure
for k=1:mp
    MB=ModalBasisP(:,k);
    for j=1:ny
        for i=1:nx
            P(i,j) = P(i,j) + sol( ndof*mx +ndof*my+ i + (k-1)*nx )*MB(j);
        end
    end
end
display(['min P =', num2str(min(min(P))),' max P =', num2str(max(max(P)))] )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=(x(2)-x(1))/2;

xvel=x(1):h:x(end);

Xv=meshgrid(xvel,Y(:,1));

X=meshgrid(x,Y(:,1));

Yvel=Y;
Yp=zeros(ny,nx);
for r=1:ne
    Yp(:,r)=Y(:,2*r-1);
end
Yp(:,nx)=Y(:,nx+ne);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( out.Contour.out )
    if(out.Contour.video)
        if(isempty(videocont))
            name='contour';
            nomefilmato=strcat(pwd,'/',name,'.avi');
            videocont= avifile(nomefilmato,'fps',5); %fps basso=filmato lento.
        end
    end
    
    if(isempty(fig_cont))
        fig_cont=figure;
        set(fig_cont,'visible',out.Contour.visible)
    else
        figure(fig_cont)
        set(fig_cont,'visible',out.Contour.visible)
    end
    clinesx=linspace(out.Cx.cinf,out.Cx.csup,20);
    clinesy=linspace(out.Cy.cinf,out.Cy.csup,20);
    clinesp=linspace(out.Cp.cinf,out.Cp.csup,20);
    
    subplot(3,1,1,'replace')
    contourf(Xv,Yvel,Ux',clinesx);
    title(['Ux t= ',num2str(t)])
    colorbar;
    ylim([-2,2])
    caxis([out.Cx.cinf,out.Cx.csup]);
    
    subplot(3,1,2,'replace')
    contourf(Xv,Yvel,Uy',clinesy);
    title(['Uy t= ',num2str(t)])
    colorbar;
    ylim([-2,2])
    caxis([out.Cy.cinf,out.Cy.csup]);
    
    subplot(3,1,3,'replace')
    contourf(X,Yp,P',clinesp);
    title(['Pressure t= ',num2str(t)])
    colorbar;
    ylim([-2,2])
    caxis([out.Cp.cinf,out.Cp.csup]);
    if(thereis(out.Contour.exp,t))
        sol=fopen(['sol_t=',num2str(t),'.out'],'w');
        Ne=(size(Xv,2)-1)/2;
        M=size(Xv,1);
        fprintf(sol,'            %u\n',3*Ne);
        fprintf(sol,'             %u\n',M);
        for k=1:Ne
            for j=1:M
                fprintf(sol,'       %1.7E       %1.7E       %1.7E\n',Xv(j,(k-1)*2 + 1),Yvel(j,(k-1)*2 + 1),Ux((k-1)*2 + 1,j));
            end
            for j=1:M
                fprintf(sol,'       %1.7E       %1.7E       %1.7E\n',Xv(j,(k-1)*2 + 2),Yvel(j,(k-1)*2 + 2),Ux((k-1)*2 + 2,j));
            end
            for j=1:M
                fprintf(sol,'       %1.7E       %1.7E       %1.7E\n',Xv(j,(k-1)*2 + 3),Yvel(j,(k-1)*2 + 3),Ux((k-1)*2 + 3,j));
            end
        end
        for k=1:Ne
            for j=1:M
                fprintf(sol,'       %1.7E       %1.7E       %1.7E\n',Xv(j,(k-1)*2 + 1),Yvel(j,(k-1)*2 + 1),Uy((k-1)*2 + 1,j));
            end
            for j=1:M
                fprintf(sol,'       %1.7E       %1.7E       %1.7E\n',Xv(j,(k-1)*2 + 2),Yvel(j,(k-1)*2 + 2),Uy((k-1)*2 + 2,j));
            end
            for j=1:M
                fprintf(sol,'       %1.7E       %1.7E       %1.7E\n',Xv(j,(k-1)*2 + 3),Yvel(j,(k-1)*2 + 3),Uy((k-1)*2 + 3,j));
            end
        end
        for k=1:Ne
            for j=1:M
                fprintf(sol,'       %1.7E       %1.7E       %1.7E\n',Xv(j,(k-1)*2 + 1),Yvel(j,(k-1)*2 + 1),P((k-1) + 1,j));
            end
            for j=1:M
                fprintf(sol,'       %1.7E       %1.7E       %1.7E\n',Xv(j,(k-1)*2 + 2),Yvel(j,(k-1)*2 + 2),(P((k-1) + 2,j)+P((k-1) + 1,j))/2);
            end
            for j=1:M
                fprintf(sol,'       %1.7E       %1.7E       %1.7E\n',Xv(j,(k-1)*2 + 3),Yvel(j,(k-1)*2 + 3),P((k-1) + 2,j));
            end
        end
    end
    
    if(out.Contour.video)
        F = export_fig(fig_cont, '-nocrop', '-a1','-transparent');
        videocont=addframe(videocont,F);
    end
    if (out.ended)
        videocont=close(videocont);
    end
end

% if( out.quiv ) %streamlines..
%     if(isempty(fig_quiver))
%         fig_quiver=figure;
%         set(fig_quiver,'visible',out.visquiver)
%     else
%         figure(fig_quiver)
%         set(fig_quiver,'visible',out.visquiver)
%     end
%     step=floor(ndof/out.nb_quiv);
%     pt=1:step:ndof;
%     quiver(Xv(:,pt),Yvel(:,pt),Ux(pt,:)',Uy(pt,:)','color','blue');
%     hold on
%     [sx,sy]=meshgrid(xvel(1),Yvel(:,1));
%     hold off
%     streamline(Xv,Yvel,Ux',Uy',sx,sy)
%     title(['Velocity field',num2str(t)])
% end

if( out.Profile.out )
    if(out.Profile.video)
        if(isempty(videoprofile))
            name='profile';
            nomefilmato=strcat(pwd,'/',name,'.avi');
            videoprofile= avifile(nomefilmato,'fps',5); %fps basso=filmato lento.
        end
    end
    if(isempty(fig_profile))
        fig_profile=figure;
        set(fig_profile,'visible',out.Profile.visible)
    else
        figure(fig_profile)
        set(fig_profile,'visible',out.Profile.visible)
    end
    ind=floor(out.Profile.x/(2*h))+1;
    plot(Yp(:,ind),P(ind,:),'linewidth',2);
    ylim([out.Cp.cinf,out.Cp.csup])
    xlim([-2,2]);
    ylabel('Pressure');
    xlabel('Y');
    title(['x= ',num2str(x(ind))])
    grid on
    if(thereis(out.Profile.exp,t))
        Pr=fopen(['P_profile_x=',num2str(x(ind)),'t=',num2str(t),'.dat'],'w');
        fprintf(Pr,'             %u\n',size(Yp,1));
        for w=1:size(Yp,1)
            fprintf(Pr,'       %1.7E       %1.7E\n',Yp(w,ind),P(ind,w));
        end
    end
    if(out.Profile.video)
        F = export_fig(fig_profile, '-nocrop', '-a1','-transparent');
        videoprofile=addframe(videoprofile,F);
    end
    if (out.ended)
        videoprofile =close(videoprofile);
    end
end

if( out.ProfileVelocity.out )
    if(out.ProfileVelocity.video)
        if(isempty(videoprofileVelocity))
            name='profileVelocity';
            nomefilmato=strcat(pwd,'/',name,'.avi');
            videoprofileVelocity= avifile(nomefilmato,'fps',5); %fps basso=filmato lento.
        end
    end
    if(isempty(fig_profileVelocity))
        fig_profileVelocity=figure;
        set(fig_profileVelocity,'visible',out.ProfileVelocity.visible)
    else
        figure(fig_profileVelocity)
        set(fig_profileVelocity,'visible',out.ProfileVelocity.visible)
    end
    ind=floor(out.ProfileVelocity.x/(h))+1;
    plot(Y(:,ind),Ux(ind,:),'linewidth',2);
    ylim([out.Cx.cinf,out.Cx.csup])
    xlim([-2,2]);
    ylabel('U_x');
    xlabel('Y');
    title(['x= ',num2str(xvel(ind))])
    grid on
    if(thereis(out.ProfileVelocity.exp,t))
        Uout=fopen(['U_profile_x=',num2str(xvel(ind)),'t=',num2str(t),'.dat'],'w');
        fprintf(Uout,'             %u\n',size(Y,1));
        for w=1:size(Y,1)
            fprintf(Uout,'       %1.7E       %1.7E\n',Y(w,ind),Ux(ind,w));
        end
    end
    if(out.ProfileVelocity.video)
        F = export_fig(fig_profileVelocity, '-nocrop', '-a1','-transparent');
        videoprofileVelocity=addframe(videoprofileVelocity,F);
    end
    if (out.ended)
        videoprofileVelocity =close(videoprofileVelocity);
    end
end

if( out.Pressure.out)
    if(out.Pressure.video)
        if(isempty(videopressure))
            name='pressure2';
            nomefilmato=strcat(pwd,'/',name,'.avi');
            videopressure= avifile(nomefilmato,'fps',5); %fps basso=filmato lento.
        end
    end
    if(isempty(fig_pressure))
        fig_pressure=figure;
        set(fig_pressure,'visible',out.Pressure.visible)
    else
        figure(fig_pressure)
        set(fig_pressure,'visible',out.Pressure.visible)
    end
    plot(X(1,:),P(:,ny/2),'linewidth',3,'color','r');
    hold on
    plot(X(1,:),P(:,3*ny/4),'linewidth',3,'color','g');
    hold on
    plot(X(1,:),P(:,ny/4),'linewidth',3,'color','k');
    hold off
    legend('center','down','up');
    ylim([out.Cp.cinf,out.Cp.csup])
    title(['Pressure at y = 0 t= ',num2str(t)])
    grid on
    
    if(thereis(out.Pressure.exp,t))
        PP=fopen(['Pressure_t=',num2str(t),'.dat'],'w');
        fprintf(PP,'             %u\n',size(X,2));
        for w=1:size(X,2)
            fprintf(PP,'       %1.7E       %1.7E\n',X(1,w),P(w,ny/2));
        end
    end
    if(out.Pressure.video)
        F = export_fig(fig_pressure, '-nocrop', '-a1','-transparent');
        videopressure=addframe(videopressure,F);
    end
    if (out.ended)
        videopressure =close(videopressure);
    end
end

if( out.Domain.out )
    if(out.Domain.video)
        if(isempty(videodomain))
            name='domain2';
            nomefilmato=strcat(pwd,'/',name,'.avi');
            videodomain= avifile(nomefilmato,'fps',5); %fps basso=filmato lento.
        end
    end
    if(isempty(fig_domain))
        fig_domain=figure;
        set(fig_domain,'visible',out.Domain.visible)
    else
        figure(fig_domain)
        set(fig_domain,'visible',out.Domain.visible)
    end
%     plot(xvel,-1+etadown,'linewidth',3);
%     hold on
%     plot(xvel,1+etaup,'linewidth',3);
%     hold off
%     ylim([-2,2])
%     title(['Computational domain at t= ',num2str(t)])
%     grid on
    
    if(thereis(out.Domain.exp,t))
        DU=fopen(['DU=',num2str(t),'.dat'],'w');
        fprintf(DU,'             %u\n',length(etadown));
        for w=1:length(etadown)
            fprintf(DU,'       %1.7E       %1.7E\n',xvel(w),-1+etadown(w));
        end
        fclose(DU);
        DD=fopen(['DD=',num2str(t),'.dat'],'w');
        fprintf(DD,'             %u\n',length(etadown));
        for w=1:length(etadown)
            fprintf(DD,'       %1.7E       %1.7E\n',xvel(w),1+etaup(w));
        end
        fclose(DD);
%         plot(xvel,-1+etadown,'linewidth',3);
%         hold on
%         plot(xvel,1+etaup,'linewidth',3);
%         hold on
%         ylim([-1.5,1.5])
%         grid on
    end
    if(out.Domain.video)
        F = export_fig(fig_domain, '-nocrop', '-a1','-transparent');
        videodomain=addframe(videodomain,F);
    end
    if (out.ended)
        videodomain =close(videodomain);
    end
end
