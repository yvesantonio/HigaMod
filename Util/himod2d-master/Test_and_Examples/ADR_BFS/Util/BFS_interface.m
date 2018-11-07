function coef=BFS_interface(uL,aRilL,bRilL,ml,mr,bc_up,bc_downL,bc_downR,Coeff_forma)
Nquad=64;
[~, yqhat, wyqhat] = gausslegendre( Nquad );
[yq,wyq]=quadrature_rule(-1,1,yqhat,wyqhat);

nxL=size(uL{1},1)/ml;

%the modal basis must be evaluated on nodes that belongs to [0,1]
MBL=new_modal_basis(ml, yq(Nquad/2+1:end), bc_up,bc_downL,Coeff_forma,1);
MBR=new_modal_basis(mr, yqhat, bc_up,bc_downR,Coeff_forma,2);

dataEvaluated=zeros(Nquad,1);
for k=1:ml
    dataEvaluated(Nquad/2+1:end)=dataEvaluated(Nquad/2+1:end)+MBL(:,k)*uL{1}( k * nxL );
end
dataEvaluated(Nquad/2+1:end)=dataEvaluated(Nquad/2+1:end)+aRilL*yq(Nquad/2+1:end)+bRilL;

coef=zeros(mr,1);
for j=1:mr
    coef(j,1)=integrate(dataEvaluated'.*MBR(:,j)',wyq);
end
% forPlot=zeros(Nquad,1);
% for j=1:mr
%     forPlot=forPlot+MBR(:,j)*coef(j,1);
% end
% %
% figure
% plot(yq,dataEvaluated)
% hold on
% plot(yq,forPlot,'r');
end