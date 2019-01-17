function [St,v]=pre_proc_VOL(St,v)

x=St.data(:,1);
y=St.data(:,2);
z=St.data(:,3);
[Tg Nm Bn]=myfrenet(x,y,z);

    [ sec j ct t]=section_mesh(St.data, v, St.tubo, 1);
    if isempty(sec)
    [ sec j ct t]=section_mesh(St.data, v, St.tubo, 2);
    end
[St.tubo v]=mesh_adapt(St.tubo,v,sec);
    [x1 y1 z1 a]=nearest(sec(1,1),sec(1,2),sec(1,3),v(:,1),v(:,2),v(:,3));
p=0;
p(1)=a;
for i=2:numel(sec)/3
    [x1 y1 z1 a]=nearest(sec(i,1),sec(i,2),sec(i,3),v(:,1),v(:,2),v(:,3));
    if a~=p(end)
        p(end+1)=a;
    end
end
p=p';
for i=1:numel(p)-1
    limit(i,:)=[p(i) p(i+1)];
end
St.lim{end+1}=limit;
clearvars limit    
St.section_start_tubo=sec;