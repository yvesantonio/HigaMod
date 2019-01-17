function [faces vertices p2]=addfaces10(faces,vertices,p1,p2);
    [vertices vp1]=addvertices(vertices,p1);
 i=1;
index=[ 1 2; 
        2 3; 
        3 1];
while i
    f=faces(i,:);
     if ismember(vp1,faces(i,:))
         for k=1:3
            t=point_in_edge(p2,vertices(f(index(k,1)),:),vertices(f(index(k,2)),:));
            if t==1
                break
            end
            
         end
         if t==1
             break
         end
     end
     i=i+1;
end
[vertices vp2]=addvertices(vertices,p2);
faces(i,:)=[vp1 vp2  f(index(k,1))];
n=numel(faces)/3;
faces(n+1,:)=[vp1 vp2  f(index(k,2))];
