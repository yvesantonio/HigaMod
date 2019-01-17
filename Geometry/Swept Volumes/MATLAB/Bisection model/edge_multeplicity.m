function multi=edge_multeplicity(f);
%compute for each edge, is multeplicity
for i=1:numel(f)/3
    k=0;
    e=[f(i,1) f(i,2)];  %e2=[f(i,1) f(i,3)];e3=[f(i,2) f(i,3)];  
     for j=1:numel(f)/3
        if (e(1)==f(j,1) | e(1)==f(j,3) | e(1)==f(j,2)) & (e(2)==f(j,1) | e(2)==f(j,3) | e(2)==f(j,2) )           
            k=k+1;
           multi(1,k,i)=j; 

        end
     end
     k=0;
     e=[f(i,2) f(i,3)];  %e2=[f(i,1) f(i,3)];e3=[f(i,2) f(i,3)];   
     for j=1:numel(f)/3 
        if (e(1)==f(j,1) | e(1)==f(j,3) | e(1)==f(j,2)) & (e(2)==f(j,1) | e(2)==f(j,3) | e(2)==f(j,2) )
            k=k+1;
           multi(2,k,i)=j;
           
        end
     end
     k=0;
     e=[f(i,1) f(i,3)];  %e2=[f(i,1) f(i,3)];e3=[f(i,2) f(i,3)];  
     for j=1:numel(f)/3
        if (e(1)==f(j,1) | e(1)==f(j,3) | e(1)==f(j,2)) & (e(2)==f(j,1) | e(2)==f(j,3) | e(2)==f(j,2) )
            k=k+1;
           multi(3,k,i)=j;
           
        end
     end
end
end