function k=spot_face(f,e)
k=0;
for j=1:numel(f)/3
    if prod(ismember(e,f(j,:)))
        k=j;
        break
    end 
end