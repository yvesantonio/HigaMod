function k=isin(point,nsec);
k=0;
for i=1:numel(nsec)/3;
    if prod(point==nsec(i,:))
        k=1;
        break
    end
end