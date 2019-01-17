function N=Nip(i,u,p,U) 
if u==1 
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