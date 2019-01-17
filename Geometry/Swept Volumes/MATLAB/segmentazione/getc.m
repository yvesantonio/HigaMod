function St=getc(C)

if numel(C)==1
St.cont=0;
St.data=C{1}(:,1:3);
St.rad=C{1}(:,4);


else
St.cont=1;
ref=C{1};
i=1;
true=1;
i_div=1;
i_eq=1;
while  true
    for j=2:numel(C)
        if C{j}(i+1,:)~=ref(i+1,:)
            div(i_div)=j;
            i_div=i_div+1;
            true=0;
        end
    end
    i=i+1;
end

Cnext1{1}=ref(i:end,:);
i_div=1;
i_eq=2;
for j=2:numel(C)
    if ismember(j,div) 
        Cnext2{i_div}=C{j}(i:end,:);
        i_div=i_div+1;
    else
        Cnext1{i_eq}=C{j}(i:end,:);
        i_eq=i_eq+1;
    end
end    
St.data=ref(1:i-1,1:3);
St.rad=ref(1:i-1,4);
St.next1=getc(Cnext1);
St.next2=getc(Cnext2);
end
end