function nf=newf(f)
k=1;

for i=1:numel(f)/3
    if prod([0 0 0]==f(i,:))
    else
        nf(k,:)=f(i,:);
        k=k+1;
    end
end
if numel(f)==k
    nf=f;
end
end
        