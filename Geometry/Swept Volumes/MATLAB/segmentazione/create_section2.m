function [St faces vertices]=create_section2(St, vertices, faces)




if St.cont==1
    St.next1=create_section(St.next1, vertices, faces);
end
if St.cont==1
    St.next2=create_section(St.next2, vertices, faces);
end
St.sec=[];


if St.cont==0
% sezione relativa alla fine del pezzo
i=numel(St.data)/3-5;
disp('inizio')
 [St.sec{end+1} n1]=section_mesh(St.data,vertices,faces,i);
St.section_end_tubo=St.sec{end};

i=1;j=1;
radn1=St.next1.rad;
radn2=St.next2.rad;
p1=St.next1.data(i,:); nel1=numel(St.next1.data)/3;
p2=St.next2.data(j,:); nel2=numel(St.next2.data)/3;

while norm(p1-p2)<1.2*(radn1(i)+radn2(j))
    if i<nel1
        i=i+1;
    end
    if j<nel2
        j=j+1;
    end
    p1=St.next1.data(i,:);
    p2=St.next2.data(j,:);
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

    
            
                
    