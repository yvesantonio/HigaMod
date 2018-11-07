function [Map,Y,W,Ynodes,Jinout]=ComputeMapFromDisplacementSecOrd(etadown_k,etaup_k,Xhat,Yold,Yolder,Yhat,FEM,DFEM,dt)

ny=size(Yhat,1); % nb of quadrature nodes in y-direction
nqnx=size(FEM,1); % nb of quadrature nodes per element in the x-direction
ne=size(Xhat(1,:),2)/nqnx; % nb of elements

B=(etadown_k+etaup_k)/2;
A=(etaup_k-etadown_k)/2;

Bquad=EvalP2(B,FEM,ne,nqnx);
DBquad=EvalP2(B,DFEM,ne,nqnx);

Aquad=EvalP2(A,FEM,ne,nqnx);
DAquad=EvalP2(A,DFEM,ne,nqnx);

Y=zeros(size(Yhat));
menoDsuJ=zeros(size(Yhat));
Ynodes=zeros(ny,2*ne+1);
for j=1:ny
    Y(j,:)=Yhat(j,1)*(1+Aquad)+Bquad;
    Ynodes(j,:)=Yhat(j,1)*(1+A)+B;
    menoDsuJ(j,:)=Yhat(j,1)*DAquad+DBquad;
end
    bordmenoDsuJ(1,:)=DAquad+DBquad;
    bordmenoDsuJ(2,:)=-DAquad+DBquad;

Jac=((1+Aquad)*ones(1,ny));
J=1./Jac;
D=-menoDsuJ'.*J;
Aup=sqrt(1+bordmenoDsuJ(1,:)'.^2);
Adown=sqrt(1+bordmenoDsuJ(2,:)'.^2);

Map=struct('Jac',Jac,'D',D,'J',J,'Aup',Aup,'Adown',Adown);


W=((3/2*Y-2*Yold+1/2*Yolder)/dt)';

Jinout=[(1+A(1))*ones(ny,1),(1+A(end))*ones(ny,1)];