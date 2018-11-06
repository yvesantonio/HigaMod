function [T N B]=myfrenet(x,y,z,init)
% Construct tangent, normal and binormal vectors of a curve
% Mehmet OZTURK - KTU Electrical and Electronics Engineering, Turkey
% "x","y",and "z" is the coordinates of the curve. The curve can be 2d
% "init" is initial normal vector that user can be specified

% % Example
% % clear,clc
% % 
% % t=2*pi*linspace(-1/2,1/2,100).';
% % 
% % x=cos(t); y=sin(t); z=t;
% % myfrenet(x,y,z)

x=x(:); y=y(:);

dim=3; sz=size(x);
if nargin==2, z=zeros(sz);  dim=2; end

z=z(:); % make column vectors

% calculate derrivatives of the curve
X=csaps(1:sz(1),x,1);
Y=csaps(1:sz(1),y,1);
Z=csaps(1:sz(1),z,1);
mx=fnval(fnder(X,1),1:sz(1)).';
my=fnval(fnder(Y,1),1:sz(1)).';
mz=fnval(fnder(Z,1),1:sz(1)).';

ind=find(sqrt(sum([mx my mz].*[mx my mz],2))>0);
data=[mx(ind) my(ind) mz(ind)]; % discard bad points

% normalize tangents
T=bsxfun(@rdivide,data,sqrt(sum(data.*data,2)));

% memory allocation for Normal and Binormal vectors
s=numel(ind);
N=zeros(s+1,3); B=zeros(s,3);

if nargin < 4 % check if user has supplied initial normal vector (init)
    % make initial normal vector perpendicular to first tangent vector
    N(1,:)=[T(1,3)*T(1,1) T(1,3)*T(1,2) -(T(1,1)^2+T(1,2)^2)];
    % check for apropirate initial normal vector (if not exist frenet normal)
    if all(N(1,:)==0) || all(isnan(N(1,:))); N(1,:)=[1 0 0]; end
else
    N(1,:)=init;
end
N(1,:)=N(1,:)./norm(N(1,:)); % normalize

% propagate the normal along the curve
for m=1:s
    B(m,:)=cross(N(m,:),T(m,:));
    N(m+1,:)=cross(T(m,:),B(m,:));
end

N(1,:)=[]; % discard initial vector
T=interp1(ind,T,1:sz(1),'*linear','extrap');
B=interp1(ind,B,1:sz(1),'*linear','extrap');
N=interp1(ind,N,1:sz(1),'*linear','extrap');

T=bsxfun(@rdivide,T,sqrt(sum(T.*T,2))); % normalize B vector
B=bsxfun(@rdivide,B,sqrt(sum(B.*B,2))); % normalize B vector
N=bsxfun(@rdivide,N,sqrt(sum(N.*N,2))); % normalize N vector

if dim==2
    T(:,3)=[]; N(:,3)=[]; B(:,3)=[];
end

if nargout==0
   figure
   if dim==3
       plot3(x,y,z,'color','r'),hold on
       quiver3(x,y,z,N(:,1),N(:,2),N(:,3),'color','g'),hold on
       quiver3(x,y,z,B(:,1),B(:,2),B(:,3),'color','b'),hold on
       grid,daspect([1 1 1]),axis vis3d
       legend('Curve','Normal','Binormal')
   else
       plot(x,y),hold on
       quiver(x,y,B(:,1),B(:,2),'color','r')
       grid,daspect([1 1 1]),axis vis3d
       legend('Tangent','Normal','Binormal')
   end
end
