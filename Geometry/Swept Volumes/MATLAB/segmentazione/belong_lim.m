function id=belong_lim(e,St)
id=0;
for j=1:numel(St.lim)
    
for i=1:numel(St.lim{j})/2
    if ismember(e,St.lim{j}(i,:))
        id=1;
        break
    end
end
if id==1
    break
end

end

if id==0 & St.cont==1
    id=belong_lim(e,St.next1);
end

if id==0 & St.cont==1
    id=belong_lim(e,St.next2);
end
    