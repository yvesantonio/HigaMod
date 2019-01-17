function [ sec i coordtg tg]=section_mesh(data, vertices, faces, i)

[Tg Nm Bn]=myfrenet(data(:,1),data(:,2),data(:,3));
 plane = [data(i,1:3), ...
          Nm(i,:), ...
          Bn(i,:)];
polyCellNow = xsecmesh(plane, vertices, faces);
while size(polyCellNow)==[1 0] | isempty(polyCellNow)
    i=i-1;
     plane = [data(i,1:3), ...
          Nm(i,:), ...
          Bn(i,:)];
polyCellNow = xsecmesh(plane, vertices, faces);
end
        min=[0 inf];
        for j=1:numel(polyCellNow)
            xm=sum(polyCellNow{j}(:,1))/numel(polyCellNow{j}(:,1));
            ym=sum(polyCellNow{j}(:,2))/numel(polyCellNow{j}(:,2));
            zm=sum(polyCellNow{j}(:,3))/numel(polyCellNow{j}(:,3));
            if norm([xm-data(i,1) ym-data(i,2) zm-data(i,3)])<min(2)
                min=[j norm([xm-data(i,1) ym-data(i,2) zm-data(i,3)])];
            end
        end
    sec=[polyCellNow{min(1)}(:,1) polyCellNow{min(1)}(:,2) polyCellNow{min(1)}(:,3)];
    disp('finitooo');
    coordtg=data(i,:);
    tg=Tg(i,:);
    