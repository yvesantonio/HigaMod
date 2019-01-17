function [St faces vertices]=create_section(St, vertices, faces)




if St.cont==1
    St.next1=create_section(St.next1, vertices, faces);
end
if St.cont==1
    St.next2=create_section(St.next2, vertices, faces);
end

St.sec=[];

if St.cont==0
    St.tuboexist=1;
    i=numel(St.data)/3-3;
    [endsec j]=section_mesh(St.data,vertices,faces,i);
    while isempty(endsec) & prod(endsec(1,:)==endsec(end,:))~=1
        i=i-1;
        [endsec j]=section_mesh(St.data,vertices,faces,i);
    end
    St.sec{end+1}=endsec;
    St.section_end_tubo=St.sec{end};
    St.fine=j;
else
    

% sezione relativa alla fine del pezzo
i=numel(St.data)/3;
disp('inizio')
[endsec j]=section_mesh(St.data,vertices,faces,i);
    while isempty(endsec) & prod(endsec(1,:)==endsec(end,:))~=1
        i=i-1;
        [endsec j]=section_mesh(St.data,vertices,faces,i);
    end
    St.sec{end+1}=endsec;
    St.section_end_tubo=St.sec{end};
    St.fine=j;
%  [s{1} n1]=section_mesh(St.next1.data,vertices,faces,i);
% [s{2} n2]=section_mesh(St.next2.data,vertices,faces,i);
% xb1=sum(s{1}(:,1))/numel(s{1}(:,1));
% yb1=sum(s{1}(:,2))/numel(s{1}(:,2));
% zb1=sum(s{1}(:,3))/numel(s{1}(:,3));
% xb2=sum(s{2}(:,1))/numel(s{2}(:,1));
% yb2=sum(s{2}(:,2))/numel(s{2}(:,2));
% zb2=sum(s{2}(:,3))/numel(s{2}(:,3));
% if norm([xb1 yb1 zb1]-St.next1.data(n1,:))<norm([xb2 yb2 zb2]-St.next2.data(n2,:))
%   St.sec{end+1}=s{1};   
%   St.fine=n1;
%   St.n=1;
% else
% St.sec{end+1}=s{2};
%   St.fine=n2;
%   St.n=2;
% end
St.section_end_tubo=St.sec{end};

i=1;j=1;
radn1=St.next1.rad;
radn2=St.next2.rad;
% vecchio metodo
% p1=St.next1.data(i,:); nel1=numel(St.next1.data)/3;
% p2=St.next2.data(j,:); nel2=numel(St.next2.data)/3;
% 
% while norm(p1-p2)<1.2*(radn1(i)+radn2(j))
%     if i<nel1
%         i=i+1;
%     end
%     if j<nel2
%         j=j+1;
%     end
%     p1=St.next1.data(i,:);
%     p2=St.next2.data(j,:);
% end
% vecchio metodo 
p11=St.next1.data(1,:); nel1=numel(St.next1.data)/3;
p22=St.next2.data(1,:); nel2=numel(St.next2.data)/3;
i=1;
while i<nel1 & norm(p11-St.next1.data(i,:))<2.5*radn1(1) 
    i=i+1;
end
j=1;
while j<nel2 & norm(p22-St.next2.data(j,:))<2.5*radn2(1) 
    j=j+1;
end

if i>=numel(St.next1.data)/3*0.85
    i=numel(St.next1.data)/3;
    St.next1.tuboexist=0;
else
    St.next1.tuboexist=1;
end

if St.next1.tuboexist
    disp('inizio')
[St.sec{end+1} i]=section_mesh(St.next1.data,vertices,faces,i);
St.next1.section_start_tubo=St.sec{end};
St.next1.start=i;
end

if j>=numel(St.next2.data)/3*0.8
    j=numel(St.next2.data)/3;
    St.next2.tuboexist=0;
else
    St.next2.tuboexist=1;
end

if St.next2.tuboexist
    disp('inizio')
[St.sec{end+1} j]=section_mesh(St.next2.data,vertices,faces,j);
St.next2.section_start_tubo=St.sec{end};
St.next2.start=j;
end

end

    
            
                
    