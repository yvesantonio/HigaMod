function [ch new_list]=chain(list)
% chain the edge to create closed profile or open profile, in case of closed
% profile last value is repeated
ch=list(1,:); % try to chain the first edge on list 
list(1,:)=[0 0];
n=numel(list)/2;
for i=1:n
    for j=1:n  
        if j~=i & ch(end)==list(j,1)
            ch(end+1)=list(j,2);
            list(j,:)=[0 0];
        end
        if j~=i & ch(end)==list(j,2)
            ch(end+1)=list(j,1);
            list(j,:)=[0 0];
        end
    end
end
k=1;
l_mod_list=n-numel(ch)+1; %first and last element coincide
for i=1:n
    if list(i,1)~=0
    new_list(k,:)=list(i,:);
    k=k+1;
    end
end
if l_mod_list==0
    new_list=-1;
else
list=new_list;
end
end
