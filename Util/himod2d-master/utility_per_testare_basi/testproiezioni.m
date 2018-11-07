clear all
close all
%**********************************************%
% Coefficienti della forma bilineare           %
%**********************************************%
chi=3;
mu    = @(x,y) (  1 + 0*x+0*y ); % Diffusione
beta1 = @(x,y) ( 20 + 0*x+0*y ); % Trasporto orizzontale
beta2 = @(x,y) (  0 + 0*x+0*y ); % Trasporto verticale
sigma = @(x,y) (  0 + 0*x+0*y ); % Reazione

Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigma,'coeffrobin',chi);

addpath ../
x=linspace(0,1,8000);
m=[2 4 8 16 32 64];% 128 256 1024  ];
dx=x(2)-x(1);
%Proiettanda
% parabola che rispetta robin
%  b=0.01;
%  mu=1;
%  C= mu*b/chi;
%  A=(-mu*b-b*chi-mu*b )/(2*mu+chi);
%  f=@(x)  A*x.^2+b*x+C;
%  df=@(x) 2*A*x+b;

% una cubica che rispetta robin
%a=-(chi+2*mu(0,0))/(3*mu(0,0)+chi);
%f=@(x) a*x.^3+x.^2;
%df=@(x) 3*a*x.^2+2*x;

% una retta
%f=@(x) x;
%df=@(x) 1+0.*x;

% un logaritmo più una retta NON CONFORME
%c=chi*log(2)*2/(mu(0,0)+2*chi);
%a=-c;
%f=@(x) c*log(1+x)+a*x;
%df=@(x) c./(1+x)+ a;
%  un esponenziale
%f=@(x) exp(x);
%df=@(x) f(x);

% un coseno
%f=@(x) cos(x);
%df=@(x) -sin(x);

%un seno periodico 1
%f=@(x) sin(2*pi*x);
%df=@(x) 2*pi*cos(2*pi*x);

%una parabola con sviluppo periodico continuo
f=@(x) x.^2-x;
df=@(x) 2*x-1;

% una cubica
% f=@(x) x.^3;
%df=@(x) 3.*x.^2;

% una parabola
% f=@(x) x.^2;
% df=@(x) 2*x;

[mb,mb_y]=new_modal_basis(m(end),x,'rob','rob',Coeff_forma);

%Con basi classiche
[mbf,mbf_y]=new_modal_basis_FOURIER(m(end),x);

%PLOT ISTRUITE
p=zeros(m(end),length(m));
p_plot=zeros(length(x),length(m));
px_plot=zeros(length(x),length(m));
for j=1:length(m)
    for i=1:m(j)
        p(i,j)      = integrami(f(x').*mb(:,i),x);
        p_plot(:,j) = mb(:,i)*p(i,j)+p_plot(:,j);
        px_plot(:,j)= mb_y(:,i)*p(i,j)+px_plot(:,j);
    end
    erroreL2(j)=sqrt(integrami( (f(x')-p_plot(:,j)).^2 ,x));
    
    erroreH1(j)=sqrt(erroreL2(j)^2+integrami(   ( df(x')-px_plot(:,j)).^2 ,x ));
end

%PLOT CLASSICHE
pf=zeros(m(end),length(m));
pf_plot=zeros(length(x),length(m));
pfx_plot=zeros(length(x),length(m));
for j=1:length(m)
    for i=1:m(j)
        pf(i,j)      = integrami(f(x').*mbf(:,i),x);
        pf_plot(:,j) = mbf(:,i)*pf(i,j)+pf_plot(:,j);
        pfx_plot(:,j)= mbf_y(:,i)*pf(i,j)+pfx_plot(:,j);
    end
    erroreL2f(j)=sqrt(integrami( (f(x')-pf_plot(:,j)).^2 ,x));
    
    erroreH1f(j)=sqrt(erroreL2f(j)^2+integrami(   ( df(x')-pfx_plot(:,j)).^2 ,x ));
end

close all
%STAMPA BASI ISTRUITE
fig=figure(1);
set(gca,...
    'FontUnits','points', ...
    'FontSize',11)
for j=1:length(m)
    plot(x,p_plot,'linewidth',2);
    hold on
end
hold on
plot(x,f(x),'k','linewidth',3)
%legend('esatta');
%title('Basi istruite')
%grid on
ylim([-0.35,0.1]);
save_img(fig,'basi_istruite','eps')
% STAMPA BASI CLASSICHE
fig=figure(2);
set(gca,...
    'FontUnits','points', ...
    'FontSize',11)
for j=1:length(m)   
    plot(x,pf_plot,'linewidth',2);
    hold on
end
hold on
plot(x,f(x),'k','linewidth',3)
%legend('esatta');
%title('Basi Classiche')
%grid on
 ylim([-0.35,0.1]);
save_img(fig,'basi_classiche','eps')
% %ERRORE ISTRUITE
% figure(3)
% loglog(m,erroreL2);
% hold on
% loglog(m,erroreH1,'r');
% 
% 
% x=log(m);
% y1=log(erroreL2);
% y2=log(erroreH1);
% 
% [a1]=polyfit(x,y1,1)
% C1=exp(a1(2));
% [a2]=polyfit(x,y2,1)
% C2=exp(a2(2));
% 
% title('Errore basi istruite')
% hold on
% loglog(m,C2*m.^-1.5,'gd')
% 
% legend('L^2','H^1','1.5')
% 
% %ERRORE CLASSICHE
% figure(4)
% loglog(m,erroreL2f);
% hold on
% loglog(m,erroreH1f,'r');
% 
% 
% x=log(m);
% y1=log(erroreL2f);
% y2=log(erroreH1f);
% 
% [a1]=polyfit(x,y1,1)
% C1=exp(a1(2));
% [a2]=polyfit(x,y2,1)
% C2=exp(a2(2));
% 
% title('Errore basi classiche')
% hold on
% loglog(m,C2*m.^-1.5,'gd')
% 
% legend('L^2','H^1','1.5')
