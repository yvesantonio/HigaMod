function isnot=not_contained(v1,v2)
isnot=1;
for i=1:numel(v2)/3
    if v1(1)==v2(i,1) & v1(2)==v2(i,2) & v1(3)==v2(i,3)
        isnot=0;
    end
end