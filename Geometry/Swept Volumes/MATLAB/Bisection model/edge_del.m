function f=edge_del(f,v,g);
i=1; n_faces=numel(f)/3;
while  i<=n_faces
    e=[f(i,1) f(i,2);
       f(i,2) f(i,3); 
       f(i,3) f(i,1)];
   for j=1:3
       u=1;
       clearvars mult
       n=numel(f)/3;
       k=1;
       while k<=n
        if ( e(j,1)==f(k,1) | e(j,1)==f(k,2) | e(j,1)==f(k,3) ) & ...
           ( e(j,2)==f(k,1) | e(j,2)==f(k,2) | e(j,2)==f(k,3) )
               multi(u)=k;
               u=u+1;
        end
        k=k+1;
       end
        if u>3
            g1=sum( v(f(multi(1),:),:) )/3; d1=norm(g1-g); g2=sum( v(f(multi(2),:),:) )/3; d2=norm(g2-g); g3=sum( v(f(multi(3),:),:) )/3; d3=norm(g3-g);
            fn=zeros(numel(f)/3-1,3);
            if max([d1,d2,d3])==d3
                d=multi(3); fn(1:d-1,:)=f(1:d-1,:); 
                fn(d:end,:)=f(d+1:end,:); 
                f=fn;
            end
            if max([d1,d2,d3])==d2 
                d=multi(2); fn(1:d-1,:)=f(1:d-1,:); fn(d:end,:)=f(d+1:end,:); f=fn;
            end
            if max([d1,d2,d3])==d1  
                d=multi(1); fn(1:d-1,:)=f(1:d-1,:); fn(d:end,:)=f(d+1:end,:); f=fn;
            end
        end
       end
    n_faces=numel(f)/3;
    i=i+1;
end
end