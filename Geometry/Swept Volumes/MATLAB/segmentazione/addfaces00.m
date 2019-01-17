function [faces, vertices]=addfaces00(faces,vertices,p1,p2)
i=1;
while i 
   %cerca la faccia che contenga il punto j-esimo e quelli j+1-esimo
   f=faces(i,:);

    v=[vertices(f(1),:); vertices(f(2),:); vertices(f(3),:)];
    index=[1 2; 2 3; 3 1];
    existp1=0;
    existp2=0;
    for k=1:3
        t=point_in_edge(p1,v(index(k,1),:),v(index(k,2),:));
        if t==1
            existp1=k;
        end
    end
    if existp1~=0
    for k=1:3
        t=point_in_edge(p2,v(index(k,1),:),v(index(k,2),:));
        if t==1
            existp2=k;
        end
    end
    end
    
if existp1~=0 & existp2~=0
    
% look for common vertex of edge
e1=[f(index(existp1,:))];
e2=[f(index(existp2,:))];

if ismember(e1(1),e2)
    cv=e1(1);  pnc1=e1(2);
else
   cv=e1(2);  pnc1=e1(1);
end

if ismember(e2(1),e1)
    cvt=e2(1);  pnc2=e2(2);
else
   cvt=e2(2);  pnc2=e2(1);
end


        break
    end

    i=i+1;   
end

[vertices ve1]=addvertices(vertices,p1);
[vertices ve2]=addvertices(vertices,p2);
faces(i,:)=[ve1 ve2 cv];
faces(end+1,:)=[ve1 ve2 pnc1];
faces(end+1,:)=[ve2 pnc2 pnc1];


