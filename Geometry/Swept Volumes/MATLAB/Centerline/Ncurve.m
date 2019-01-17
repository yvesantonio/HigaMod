function C=Ncurve(Px,Py,Pz,el,deg)
% el numero punti della nuova curca, equidistanziati
% if 2 dimension curve set the coordinate whic is not used to costant value
% deg is degree of spline
num=numel(Px);
P=[Px,Py,Pz];
U=zeros(1,num+deg+1);
U(end-deg:end)=1;
mid=numel(U)-2*deg-2;
for i=1:mid
    U(i+deg+1)=1/(mid+1)*i;
end
clearvars Cx Cy Cz
u=linspace(0,0.99999,el);
for j=1:numel(u)
    Cx(j)=0;
    Cy(j)=0;
    Cz(j)=0;
        for i=1:num
            Cx(j)=Cx(j)+Nip(i,u(j),deg,U)*P(i,1);    
            Cy(j)=Cy(j)+Nip(i,u(j),deg,U)*P(i,2);  
            Cz(j)=Cz(j)+Nip(i,u(j),deg,U)*P(i,3);
        end
end
C=[Cx;Cy;Cz]';
end

function N=Nip(i,u,p,U) 
if u==U(end)
    N=1;
else

    if p==0
        if u>=U(i) & u<U(i+1)
            N=1;
        else
            N=0;
        end
    else
        if U(i+p)-U(i) == 0
            N1=0;
        else
            N1=(u-U(i))/(U(i+p)-U(i))*Nip(i,u,p-1,U);
        end
        if U(i+p+1)-U(i+1) == 0
            N2=0;
        else
            N2=(U(i+p+1)-u)/(U(i+p+1)-U(i+1))*Nip(i+1,u,p-1,U);
        end
        N=N1+N2;
    end 
end
end