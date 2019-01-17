%% sorting of the senterline to get separated curve
mystlPlot(vertices,faces,0.3);hold on
plot3(x,y,z)
mystlPlot(vertices,faces,0.3);hold on
for i=1:numel(x)
plot3(x(i),y(i),z(i),'d','Color',[1-i/numel(x) 0 i/numel(x)])
pause (0.05)
end

%% separete the 2 centerline
clearvars C1 C2 C3

C1(1,:)=[x(1) y(1) z(1)];
jump=norm([x(1) y(1) z(1)]-[x(2) y(2) z(2)]);
k=2;
while jump<radius_max(k-1)
    C1(k,:)=[x(k) y(k) z(k)];
    jump=norm([x(k) y(k) z(k)]-[x(k+1) y(k+1) z(k+1)]);
    k=k+1;
end

not_in=1;
k=k+1;
i=2;
C2(1,:)=[x(k) y(k) z(k)];
not_in=not_contained([x(k+1) y(k+1) z(k+1)],C1);
while not_in
    C2(i,:)=[x(k+1) y(k+1) z(k+1)];
    i=i+1;
    k=k+1;
    not_in=not_contained([x(k+1) y(k+1) z(k+1)],[C1]);
end
C3=[x(k:end) y(k:end) z(k:end)];

% C1 modification
not_in=1;
i=1;
while not_in
    not_in=not_contained([C1(i,:)],[C3]);
    i=i+1;
end
C1=C1(1:i-1,:);
 
mystlPlot(vertices,faces,0.3);hold on
% plot3(x,y,z,'b')
plot3(C1(:,1),C1(:,2),C1(:,3),'g',C2(:,1),C2(:,2),C2(:,3),'r',C3(:,1),C3(:,2),C3(:,3),'b')