function [f v]=mesh_adapt(f,v,sec)

index=[1 2;2 3; 3 1];
for n=1:numel(sec)/3
    p=sec(n,:);
    dist(n,1)=inf;
for i=1:numel(f)/3
    fi=f(i,:);
    for ii=1:3
        v1=v(fi(index(ii,1)),:); v2=v(fi(index(ii,2)),:);
        d=norm(v2-v1); pm=(v2+v1)/2;
         if norm(pm-p)<d
            dist_new=norm(cross(p-v1,v2-v1))/norm(v2-v1);
            if dist_new<dist(n,1)
                dist(n,1)=dist_new;
                edge(n,:)=fi(index(ii,:));
                plot3([v1(1) v2(1)],[v1(2) v2(2)] ,[v1(3) v2(3)] ,'r')
             end
        
        end
    end
end
end
faccia=zeros(size(edge));
for i=1:numel(sec)/3
    k=1;
    for ii=1:numel(f)/3
        if ismember(edge(i,:),f(ii,:))
            faccia(i,k)=ii;
            k=k+1;
            if k==3
                break
            end
        end
    end
end
%% inserimento mesh
for i=1:numel(sec)/3-1
    
    if ismember(faccia(i,1),faccia(i+1,:))
        cf=faccia(i,1);
    else
        cf=faccia(i,2);
    end
    
    if ismember(edge(i,1),edge(i+1,:))
        cv=edge(i,1); ncv1=edge(i,2);
    else
        cv=edge(i,2); ncv1=edge(i,1);
    end
    
    if ismember(edge(i+1,1),edge(i,:))
        ncv2=edge(i+1,2);
    else
        ncv2=edge(i+1,1);
    end
    [v, vp1]=addvertices(v,sec(i,:));
    [v, vp2]=addvertices(v,sec(i+1,:));
    if cf~=0
    f(cf,:)=[cv vp1 vp2];
    f(end+1,:)=[vp1 vp2 ncv1];
    f(end+1,:)=[ncv1 ncv2 vp2];
    end
    
    
    
    
    
    
    
end