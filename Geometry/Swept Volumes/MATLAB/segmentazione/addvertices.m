function [v add]=addvertices(v,p)
add=0;
for i=1:numel(v)/3
    if v(i,:)==p
        add=i;
        break
    end
    i=i+1;
end
if add==0
v(end+1,:)=p;
add=numel(v)/3;
end
end
        