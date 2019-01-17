function [nf f]=segmento(nf,f,e,i,St,v)
f_temp=f(i,:);
f(i,:)=[0 0 0];
A=[f_temp(1) f_temp(2);
   f_temp(2) f_temp(3); 
   f_temp(3) f_temp(1)];
e1=[e(1) 0];
e2=[e(2) 0];
for k=1:3
    if A(k,1)~=e(1) &  A(k,1)~=e(2)
        e1(2)=A(k,1);
        e2(2)=A(k,1);
        break
    end
end
j1=spot_face(f,e1);

lim=belong_lim(e1,St);

if j1~=0 &  ~lim
    [nf f]=segmento(nf,f,e1,j1,St,v);
end

lim=belong_lim(e2,St);

j2=spot_face(f,e2);
if j2~=0 &  ~lim
    [nf f]=segmento(nf,f,e2,j2,St,v);
end
disp('a')
nf(end+1,:)=A(:,1)';


end
