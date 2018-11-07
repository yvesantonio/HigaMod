function sol=evalHiMod(u,mb,nodex,numquadnode,hx)
sol=zeros(length(nodex),size(mb,1));
ne=length(nodex)/numquadnode;
nx=ne+1;

for h=1:ne
    xl=(h-1)*hx;
    xr=h*hx;
    x1=nodex((h-1)*numquadnode +1);
    x2=nodex((h-1)*numquadnode +2);
    x3=nodex((h-1)*numquadnode +3);
    x4=nodex((h-1)*numquadnode +4);
    for k=1:size(mb,1)
        for imb=1:size(mb,2)
            
            ul=u(h+(imb-1)*nx);
            ur=u(h+1+(imb-1)*nx);
                       
            c1=1/hx*(ur*(x1-xl)-ul*(x1-xr));
            sol((h-1)*numquadnode +1,k)=sol((h-1)*numquadnode +1,k) + c1*mb(k,imb);
            
            c2=1/hx*(ur*(x2-xl)-ul*(x2-xr));
            sol((h-1)*numquadnode +2,k)=sol((h-1)*numquadnode +2,k) + c2*mb(k,imb);
            
            c3=1/hx*(ur*(x3-xl)-ul*(x3-xr));
            sol((h-1)*numquadnode +3,k)=sol((h-1)*numquadnode +3,k) + c3*mb(k,imb);
            
            c4=1/hx*(ur*(x4-xl)-ul*(x4-xr));
            sol((h-1)*numquadnode +4,k)=sol((h-1)*numquadnode +4,k) + c4*mb(k,imb);
        end
    end
end