function coef=BFS_interfaceNeu(uR,ml,mr,bc_up,bc_downL,bc_downR,Coeff_forma,hx)
Nquad=64;
[~, yq, wyq] = gausslegendre( Nquad );

nxR=size(uR{1},1)/mr;

%the modal basis must be evaluated on nodes that belongs to [0,1]
MBL=new_modal_basis(ml, yq, bc_up,bc_downL,Coeff_forma,1);
MBR=new_modal_basis(mr, yq/2+0.5, bc_up,bc_downR,Coeff_forma,2);

dataEvaluated1=zeros(Nquad,1);
for k=1:mr
    dataEvaluated1=dataEvaluated1+MBR(:,k)*uR{1}( (k-1) * nxR + 1 );
end
dataEvaluated2=zeros(Nquad,1);
for k=1:mr
    dataEvaluated2=dataEvaluated2+MBR(:,k)*uR{1}( (k-1) * nxR + 2 );
end
dataEvaluated=(dataEvaluated2-dataEvaluated1)/hx;

coef=zeros(ml,1);
for j=1:ml
    coef(j,1)=integrate(dataEvaluated'.*MBL(:,j)',wyq);
end
% forPlot=zeros(Nquad,1);
% for j=1:ml
%     forPlot=forPlot+MBL(:,j)*coef(j,1);
% end
% 
% figure
% plot(yq,dataEvaluated)
% hold on
% plot(yq,forPlot,'r');
end