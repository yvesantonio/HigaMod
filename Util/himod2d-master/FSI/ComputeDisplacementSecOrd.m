function [etadown_k,etaup_k]=ComputeDisplacementSecOrd(u,dt,etadown_old,etaup_old,etadown_older,etaup_older,MB)
%
% It computes the new displacement given the velocity in y-direction and
% the displacement from the previous time step
%
fedof=length(etadown_old);
etaup_k=zeros(fedof,1);
etadown_k=etaup_k;
for k=1:size(MB,2)
        etaup_k   = etaup_k    + 2/3*dt*u( (k-1)*fedof + 1 : k*fedof )*MB(1,k);
        etadown_k = etadown_k  + 2/3*dt*u( (k-1)*fedof + 1 : k*fedof )*MB(2,k);
end
etaup_k=etaup_k+4/3*etaup_old-1/3*etaup_older;
etadown_k=etadown_k+4/3*etadown_old-1/3*etadown_older;