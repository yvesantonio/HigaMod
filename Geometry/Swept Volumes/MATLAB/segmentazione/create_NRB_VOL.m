function St=create_NRB_VOL(St,v)
if St.cont==1
    St.next1=create_NRB_VOL(St.next1,v);
    St.next2=create_NRB_VOL(St.next2,v);
end
if St.tuboexist==1
disp('fffffiga')
[faces v]=closefaces(St.tubo,v);
fine=St.fine;
start=St.start;

x=St.data(:,1);
y=St.data(:,2);
z=St.data(:,3);
[Tg Nm Bn]=myfrenet(x,y,z);
for i=1:numel(x)-1
    if Tg(i,:)*Tg(i+1,:)'<0
        Tg(i+1,:)=-Tg(i+1,:);
    end
end
% if start~=1
srf_pti{1}=St.section_start_tubo;
plot3(srf_pti{1}(:,1),srf_pti{1}(:,2),srf_pti{1}(:,3),'k','LineWidth',2)       
Ctg(1,:)=[x(start) y(start) z(start)];
tg(1,:)=Tg(start,:);
k=2;

span=start+4:4:St.fine-4;
if span(end)~=St.fine-4
    span(end+1)=St.fine-4;
end
    
for i=span
    j=i;
    [ srf_pti{k} j ct t]=section_mesh(St.data, v, faces, j);
    plot3(srf_pti{k}(:,1),srf_pti{k}(:,2),srf_pti{k}(:,3),'k','LineWidth',2)
    index(i)=j;
    Ctg(k,:)=ct;
    tg(k,:)=Tg(j,:);
    k=k+1;
end
    srf_pti{end+1}=St.section_end_tubo;
Ctg(end+1,:)=[x(St.fine) y(St.fine) z(St.fine)];
tg(end+1,:)=Tg(St.fine,:);
plot3(srf_pti{end}(:,1),srf_pti{end}(:,2),srf_pti{end}(:,3),'k','LineWidth',2)


 srf_pti=rotate_and_ordinate(srf_pti,Ctg,tg);

    
%% get coons surface in nurbs structure
for i=1:numel(srf_pti)
    coons{i}=get_coons(srf_pti{i},6,2);
end
%%
k=1;
for i = 1:numel(srf_pti)
    surf_i = coons{i};
    volCoefs(:,:,:,k) = surf_i.coefs;
    k=k+1;
end
%% deg of volume 
d=2;
%% on U direction
n=numel(volCoefs(1,:,1,1));
U=zeros(1,n+d+1); % automatic node creation
U(end-d:end)=1;
mid=numel(U)-2*d-2;
for i=1:mid
    U(i+d+1)=1/(mid+1)*i;
end
%% on V direction
n=numel(volCoefs(1,1,:,1));
V=zeros(1,n+d+1); % automatic node creation
V(end-d:end)=1;
mid=numel(V)-2*d-2;
for i=1:mid
    V(i+d+1)=1/(mid+1)*i;
end
%% on W direction
n=numel(volCoefs(1,1,1,:));
W=zeros(1,n+d+1); % automatic node creation
W(end-d:end)=1;
mid=numel(W)-2*d-2;
for i=1:mid
    W(i+d+1)=1/(mid+1)*i;
end
%%
St.VOL = nrbmak(volCoefs,{U,V,W});
% figure
% nrbplot(St.VOL,[10 10 50]);
disp('fine creazione di un VOL')
% pause()
end

    