function d=DerNip(k,i,u,U,p)
if p==0
    d=0;
else
    if k==1 %fist derivative
        if U(i+p)-U(i)==0
        N1=0;
        else
        N1=p*Nip(i,u,p-1,U)/(U(i+p)-U(i));     
    end
    if U(i+p+1)-U(i+1)==0
        N2=0;
        else
        N2=p*Nip(i+1,u,p-1,U)/(U(i+p+1)-U(i+1));
    end
    
    d=N1-N2;
    else % derivative k-th, k>1
        
    if U(i+p)-U(i)==0
        N1=0;
        else
        N1=DerNip(k-1,i,u,U,p-1)/(U(i+p)-U(i));
    end
    if U(i+p+1)-U(i+1)==0
        N2=0;
        else
        N2=DerNip(k-1,i+1,u,U,p-1)/(U(i+p+1)-U(i+1));
    end    
    d=p*(N1-N2);
    end
end
end

    