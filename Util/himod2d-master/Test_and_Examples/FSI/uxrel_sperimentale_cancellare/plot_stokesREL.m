function   plot_stokes(sol,x,Y,ModalBasisX,ModalBasisY,ModalBasisP,out,t,etadown,etaup,mx,my,mp)

persistent fig_cont fig_quiver fig_profile fig_pressure fig_domain
persistent videopressure videocont videoprofile videodomain

nx=length(x);
ne=nx-1;
ndof=ne+nx;
ny=length(Y(:,1));
Ux=zeros(ndof,ny);
Uy=Ux;
P=zeros(nx,ny);

%X- velocity
for k=3:mx+2
    MB=ModalBasisX(:,k);
    for j=1:ny
        for i=1:ndof
            Ux(i,j) = Ux(i,j) + sol( i + (k-1-2)*ndof+2 )*MB(j);
        end
    end
end
for k=1:2
    MB=ModalBasisX(:,k);
    for j=1:ny
            Ux(1,j) = Ux(1,j) + sol( k )*MB(j);
    end
end
display(['min Ux =', num2str(min(min(Ux))),' max Ux =', num2str(max(max(Ux)))] )

%Y- velocity
for k=1:my
    MB=ModalBasisY(:,k);
    for j=1:ny
        for i=1:ndof
            Uy(i,j) = Uy(i,j) + sol(2 +ndof*mx + i + (k-1)*ndof )*MB(j);
        end
    end
end
display(['min Uy =', num2str(min(min(Uy))),' max Uy =', num2str(max(max(Uy)))] )


%Pressure
for k=1:mp
    MB=ModalBasisP(:,k);
    for j=1:ny
        for i=1:nx
            P(i,j) = P(i,j) + sol( 2+ndof*mx +ndof*my+ i + (k-1)*nx )*MB(j);
        end
    end
end
display(['min P =', num2str(min(min(P))),' max P =', num2str(max(max(P)))] )

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

if( out.cont )
    if(out.videocont)
        if(isempty(videocont))
            name='contour';
            nomefilmato=strcat(pwd,'/',name,'.avi');
            videocont= avifile(nomefilmato,'fps',5); %fps basso=filmato lento.
        end
    end
    if(isempty(fig_cont))
        fig_cont=figure;
        set(fig_cont,'visible',out.viscont)
    else
        figure(fig_cont)
        set(fig_cont,'visible',out.viscont)
    end
    
    
    clinesx=linspace(out.Cx.cinf,out.Cx.csup,20);
    clinesy=linspace(out.Cy.cinf,out.Cy.csup,20);
    clinesp=linspace(out.Cp.cinf,out.Cp.csup,20);
    
    plot(Y(:,1),Ux(1,:),'linewidth',2);
    % for the womersley test
    %     hold on
    %     plot(Y(:,1),exactsol(Y(:,1),t),'linewidth',3,'color','red')
    %     hold off
%     ylim([out.Cx.cinf,out.Cx.csup])
%     xlim([-2,2]);
%     title(['Ux at x=',num2str(xvel(1)),' t= ',num2str(t)])
%     ylabel('Y');
    grid on
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
        
    if(out.videocont)
        F = export_fig(fig_cont, '-nocrop', '-a1','-transparent');
        videocont=addframe(videocont,F);
    end
    if (out.ended)
        videocont=close(videocont);
    end
end

if( out.quiv )
    if(isempty(fig_quiver))
        fig_quiver=figure;
        set(fig_quiver,'visible',out.visquiver)
    else
        figure(fig_quiver)
        set(fig_quiver,'visible',out.visquiver)
    end
    step=floor(ndof/out.nb_quiv);
    pt=1:step:ndof;
    quiver(Xv(:,pt),Yvel(:,pt),Ux(pt,:)',Uy(pt,:)','color','blue');
    hold on
    [sx,sy]=meshgrid(xvel(1),Yvel(:,1));
    hold off
    streamline(Xv,Yvel,Ux',Uy',sx,sy)
    title(['Velocity field',num2str(t)])
end
if( out.profile )
    if(out.videoprofile)
        if(isempty(videoprofile))
            name='profile';
            nomefilmato=strcat(pwd,'/',name,'.avi');
            videoprofile= avifile(nomefilmato,'fps',5); %fps basso=filmato lento.
        end
    end
    if(isempty(fig_profile))
        fig_profile=figure;
        set(fig_profile,'visible',out.visprofile)
    else
        figure(fig_profile)
        set(fig_profile,'visible',out.visprofile)
    end
    
    plot(Y(:,1),P(1,:),'linewidth',2);
    % for the womersley test
    %     hold on
    %     plot(Y(:,1),exactsol(Y(:,1),t),'linewidth',3,'color','red')
    %     hold off
    ylim([out.Cp.cinf,out.Cp.csup])
    xlim([-2,2]);
    %title(['Ux at x=',num2str(xvel(floor(ndof/2))),' t= ',num2str(t)])
    ylabel('Y');
    grid on
    %tol=1e-10;
    %     if(abs(t-1)<tol||abs(t-1.4)<tol||abs(t-1.5)<tol)%||abs(t-1.6)<tol||abs(t-2)<tol||abs(t-2.3)<tol)
    %         export_fig(fig_profile, 'filename',['t=',num2str(t)],'-png','-transparent')
    %     end
    if(out.videoprofile)
        F = export_fig(fig_profile, '-nocrop', '-a1','-transparent');
        videoprofile=addframe(videoprofile,F);
    end
    if (out.ended)
        videoprofile =close(videoprofile);
    end
end
if( out.pressure )
    if(out.videopressure)
        if(isempty(videopressure))
            name='pressure';
            nomefilmato=strcat(pwd,'/',name,'.avi');
            videopressure= avifile(nomefilmato,'fps',5); %fps basso=filmato lento.
        end
    end
    if(isempty(fig_pressure))
        fig_pressure=figure;
        set(fig_pressure,'visible',out.vispressure)
    else
        figure(fig_pressure)
        set(fig_pressure,'visible',out.vispressure)
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
    
    if(out.videopressure)
        F = export_fig(fig_pressure, '-nocrop', '-a1','-transparent');
        videopressure=addframe(videopressure,F);
    end
    if (out.ended)
        videopressure =close(videopressure);
    end
end

if( out.domain )
    if(out.videodomain)
        if(isempty(videodomain))
            name='domain';
            nomefilmato=strcat(pwd,'/',name,'.avi');
            videodomain= avifile(nomefilmato,'fps',5); %fps basso=filmato lento.
        end
    end
    if(isempty(fig_domain))
        fig_domain=figure;
        set(fig_domain,'visible',out.visdomain)
    else
        figure(fig_domain)
        set(fig_domain,'visible',out.visdomain)
    end
    plot(xvel,-1+etadown,'linewidth',3);
    hold on
    plot(xvel,1+etaup,'linewidth',3);
    hold off
    ylim([-2,2])
    title(['Computational domain at t= ',num2str(t)])
    grid on
    
    if(out.videodomain)
        F = export_fig(fig_domain, '-nocrop', '-a1','-transparent');
        videodomain=addframe(videodomain,F);
    end
    if (out.ended)
        videodomain =close(videodomain);
    end
end