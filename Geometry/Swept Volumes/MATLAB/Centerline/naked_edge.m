function n_e=naked_edge(f);
% f is the face vector f(i,:)= Vx Vy Vz
% list is the list of the nakerd edge
% list(i,:)= V1 V2
% get back -1 -1 is no naked edge exist :)
n_e(1,:)=[-1 -1];
k=1;
n=numel(f)/3;
for i=1:n
    e=[f(i,1) f(i,2)];  %e2=[f(i,1) f(i,3)];e3=[f(i,2) f(i,3)];  
    isnaked=1;  
    j=1;
     while j<=n & isnaked
        if j~=i & (e(1)==f(j,1) | e(1)==f(j,3) | e(1)==f(j,2)) & (e(2)==f(j,1) | e(2)==f(j,3) | e(2)==f(j,2) )
           
       isnaked=0;
        end
       j=j+1;
     end
     if isnaked==1
         n_e(k,:)=e;
         k=k+1;
     end 
     
     e=[f(i,2) f(i,3)];  %e2=[f(i,1) f(i,3)];e3=[f(i,2) f(i,3)];  
    isnaked=1;  
    
    j=1;
     while j<=n & isnaked
        if j~=i & (e(1)==f(j,1) | e(1)==f(j,3) | e(1)==f(j,2)) & (e(2)==f(j,1) | e(2)==f(j,3) | e(2)==f(j,2) )
           
       isnaked=0;
        end
       j=j+1;
     end
     if isnaked==1
         n_e(k,:)=e;
         k=k+1;
     end 
     
     e=[f(i,1) f(i,3)];  %e2=[f(i,1) f(i,3)];e3=[f(i,2) f(i,3)];  
    isnaked=1;  
    j=1;
     while j<=n & isnaked
        if j~=i & (e(1)==f(j,1) | e(1)==f(j,3) | e(1)==f(j,2)) & (e(2)==f(j,1) | e(2)==f(j,3) | e(2)==f(j,2) )
           
       isnaked=0;
        end
       j=j+1;
     end
     if isnaked==1
         n_e(k,:)=e;
         k=k+1;
     end
end

end