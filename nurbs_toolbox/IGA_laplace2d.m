%-+--------------------------------------------------------------------
% HiModlab is a general purpose Hierarchical model reduction library.
% Copyright (C) 2006-2017  by the HiMod authors (see authors.txt).
%
% This file is part of the HiModlab library.
%
% The HiModlab library is free software: you can use it, redistribute
% it and/or modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation, either
% version 2 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Acknowledgement and Disclaimer
%
% The development of HiModLab was supported by the US National Science
% Foundation (NSF Project, DMS 1419060 ). Any opinions, findings and conclusions
% or recomendations expressed in the source code and documentation are
% those of the authors and do not necessarily reflect the views of the
% National Science Foundation (NSF).
%
%
%
%-+--------------------------------------------------------------------
function error = IGA_laplace2d(p,q,mcp,ncp)

% A (very simple) 2D IGA code
%
% Solve -L(u) = f, where L is the Laplace operator, with homogeneous Dirichlet b.c.'s
%        and with a manufactured body load s.t. the analytical solution is available
% Domain: quarter of an annulus
% 
% Input:  p,q     = approximation degrees in the two parametric directions
%         mcp,ncp = number of control points in the two parametric directions
%
% Output: error   = (approximation of the) L2-norm error
%
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % By:                                            %
%         %          Alessandro Reali                      %
%         %          University of Pavia                   %
%         % email:   alessandro.reali@unipv.it             %
%         % webpage: http://www.unipv.it/alereali          %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set number of Gauss points for integration
ngaussi = p + 1;
if p > 9
    ngaussi = 15;
    display('ngaussi = 15... check Gauss rule for very high degrees!')
end

ngaussj = q + 1;
if q > 9
    ngaussj = 15;
    display('ngaussj = 15... check Gauss rule for very high degrees!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial geometry: quarter annulus (internal & external radii: Ri & Re)
% [to change geometry just insert below the initial control net and knot vectors]

% set initial control net (i.e., control point coordinates and weights)
Ri = 1;
Re = 3;
X = [Ri Re;
     Ri Re;
     0  0];
Y = [0  0;
     Ri Re;
     Ri Re];
[mcp0,ncp0] = size(X);
w = ones(mcp0,ncp0);
w(2,:) = sqrt(2)/2;

% set initial degrees and knot vectors
p0 = 2;
q0 = 1;
u_knot = augknt([0 1],p0+1);
v_knot = augknt([0 1],q0+1);

BB = zeros(4,mcp0,ncp0);
BB(1,:,:) = X.*w;
BB(2,:,:) = Y.*w;
BB(4,:,:) = w;
nurbs = nrbmak(BB,{u_knot v_knot});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p-refinement
inc_p = p - p0;
inc_q = q - q0;
nurbs = nrbdegelev(nurbs,[inc_p inc_q]);
% h-refinement
nspani = mcp - p;
nspanj = ncp - q;
hi = 1/nspani;
hj = 1/nspanj;
nurbs = nrbkntins(nurbs,{hi:hi:1-hi, hj:hj:1-hj});

u_knot = nurbs.knots{1};
v_knot = nurbs.knots{2};

BB = nurbs.coefs;
Bx = squeeze(BB(1,:,:)./BB(4,:,:));
By = squeeze(BB(2,:,:)./BB(4,:,:));
w  = squeeze(BB(4,:,:));

% intializations
nel = nspani*nspanj;
nnod = mcp*ncp;
nsh = (p+1)*(q+1);

f_gl = zeros(nnod,1);
row = zeros(1,nnod^2); 
col = zeros(1,nnod^2); 
val = zeros(1,nnod^2);
icount = 0;

% set element connectivity
[ien,inn] = gen_ien_inn(nel,p,q,mcp,ncp,nnod,nsh);

% set data for homogeneous Dirichlet b.c.'s
bc = ones(mcp,ncp);
bc(2:mcp-1,2:ncp-1) = 0;
bc = reshape(bc,1,nnod);
ibc = find(bc);
clear bc

% precompute values of univariate shape functions and derivatives at Gauss points
[gpi,gwi] = gauss(ngaussi);
[gpj,gwj] = gauss(ngaussj);

gaussi = zeros(1,ngaussi*nspani);
gaussj = zeros(1,ngaussj*nspanj);
for ni = p+1:p+nspani;
    gaussi((ni-p-1)*ngaussi+1:(ni-p)*ngaussi) =... 
           ((u_knot(ni+1) - u_knot(ni))*gpi + u_knot(ni+1) + u_knot(ni))/2;
end
for nj = q+1:q+nspanj;
    gaussj((nj-q-1)*ngaussj+1:(nj-q)*ngaussj) =... 
           ((v_knot(nj+1) - v_knot(nj))*gpj + v_knot(nj+1) + v_knot(nj))/2;
end
M = spcol(u_knot,p+1,sort([gaussi,gaussi]));
N = spcol(v_knot,q+1,sort([gaussj,gaussj]));
clear gaussi gaussj

% loop over elements for local array construction
for iel = 1:nel
    ni = inn(ien(iel,1),1);
    nj = inn(ien(iel,1),2);
    if((u_knot(ni) ~= u_knot(ni+1)) && (v_knot(nj) ~= v_knot(nj+1)))
        da = (u_knot(ni+1) - u_knot(ni))*(v_knot(nj+1) - v_knot(nj))/4;
        k_el = zeros(nsh,nsh);
        f_el = zeros(nsh,1);
        % loop over gauss points
        for igauss = 1:ngaussi
            for jgauss = 1:ngaussj
                Mi = M(2*((ni-p-1)*ngaussi+igauss)-[1,0],ni:-1:ni-p);
                Nj = N(2*((nj-q-1)*ngaussj+jgauss)-[1,0],nj:-1:nj-q);
                Bxij = reshape(Bx(ni:-1:ni-p,nj:-1:nj-q),nsh,1);
                Byij = reshape(By(ni:-1:ni-p,nj:-1:nj-q),nsh,1);
                wij  = reshape(w(ni:-1:ni-p,nj:-1:nj-q),nsh,1);
                
                [shape,dshape,J,xx,yy] = eval_sh_fcn(Mi,Nj,Bxij,Byij,wij,nsh);
                gwt = gwi(igauss)*gwj(jgauss)*J*da;
                
                % manufactured body load
                ff = (2*xx^4 - 50*xx^2 - 50*yy^2 + 2*yy^4 + 4*xx^2*yy^2 + 100)*sin(xx)*sin(yy) +... 
                     (68*xx - 8*xx^3 - 8*xx*yy^2)*cos(xx)*sin(yy) +... 
                     (68*yy - 8*yy^3 - 8*yy*xx^2)*cos(yy)*sin(xx);
                
                % compute local r.h.s. vector and stiffness 
                f_el = f_el + shape*ff*gwt;
                k_el = k_el + (dshape(:,1)*dshape(:,1)' + dshape(:,2)*dshape(:,2)')*gwt;
            end
        end
        
        % assemble r.h.s. vector and prepare vectors for sparse stiffness assembly
        for a = 1:nsh
            i = ien(iel,a);
            if (i ~= ibc)
                f_gl(i) = f_gl(i) + f_el(a);
                for b = 1:nsh
                    j = ien(iel,b);
                    icount = icount + 1;
                    row(icount) = i;
                    col(icount) = j;
                    val(icount) = k_el(a,b);
                end
            else % apply homogeneous Dirichlet b.c.'s 
                icount = icount + 1;
                row(icount) = i;
                col(icount) = i;
                val(icount) = 1;
            end
        end
    end
end
clear inn ien

row = row(1:icount);
col = col(1:icount);
val = val(1:icount);

% assemble stiffness
k_gl = sparse(row,col,val,nnod,nnod);
clear row col val

% solve linear system
sol = k_gl\f_gl;
clear k_gl f_gl


% evaluate relative L^2-norm error versus the manufactured analytical solution
BB(3,:,:) = reshape(sol,mcp,ncp).*w;
nurbs = nrbmak(BB,{u_knot v_knot});

pnts = nrbeval(nurbs,{linspace(0,1,max(100,mcp)) linspace(0,1,max(100,ncp))});
xx = squeeze(pnts(1,:,:));
yy = squeeze(pnts(2,:,:));
approx = squeeze(pnts(3,:,:));
exact = (xx.^2 + yy.^2 - 1).*(xx.^2 + yy.^2 - 16).*sin(xx).*sin(yy);
error = sqrt(sum(sum((exact-approx).^2))/sum(sum(exact.^2)));

% plot approximate solution
figure(1)
surf(xx,yy,approx,'facecolor','interp','edgecolor','none','facelighting','phong')
title('approximate solution') 
view(0,90), colorbar, axis square

% plot difference between exact and approximate solution
figure(2)
surf(xx,yy,exact-approx,'facecolor','interp','edgecolor','none','facelighting','phong')
title('difference between exact and approximate solution') 
view(0,90), colorbar, axis square

return

    
function [ien,inn] = gen_ien_inn(nel,p,q,mcp,ncp,nnod,nsh)

inn = zeros(nnod,2);
ien = zeros(nel,nsh);
g = 0;
e = 0;
for j = 1:ncp
    for i = 1:mcp
        g = g + 1;
        inn(g,1) = i;
        inn(g,2) = j;
        if((i>=p+1) && (j>=q+1))
            e = e + 1;
            for k = 0:q
                for l = 0:p
                    ien(e,(p+1)*k+l+1) = g - mcp*k - l;
                end
            end
        end
    end
end  

return    


function [shape,dshape,J,xx,yy] = eval_sh_fcn(M,N,Bx,By,w,nsh)

shape  = reshape((M(1,:)'*N(1,:)),nsh,1).*w;
dshapex = reshape((M(2,:)'*N(1,:)),nsh,1).*w;
dshapey = reshape((M(1,:)'*N(2,:)),nsh,1).*w;

denom_sum  = sum(shape);
derv_sum_u = sum(dshapex);
derv_sum_v = sum(dshapey);

dshapex = dshapex/denom_sum - (shape*derv_sum_u)/(denom_sum^2);
dshapey = dshapey/denom_sum - (shape*derv_sum_v)/(denom_sum^2);
shape   = shape/denom_sum;

xx = shape'*Bx;
yy = shape'*By;

dxdxi      = zeros(2,2);
dxdxi(1,1) = dshapex'*Bx;
dxdxi(1,2) = dshapey'*Bx;
dxdxi(2,1) = dshapex'*By;
dxdxi(2,2) = dshapey'*By;

J = dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1);

dshape = [dshapex,dshapey]/dxdxi;

return


function [GP,GW] = gauss(ng)

if (ng==1)
   GP(1) = 0;
   GW(1) = 2;
elseif (ng==2)
   GP(1) = -.5773502691896259;
   GP(2) = .5773502691896259;
   GW(1) = 1;
   GW(2) = 1;
elseif (ng==3)
   GP(1) = -.7745966692414834;
   GP(2) = 0;
   GP(3) = .7745966692414834;
   GW(1) = .5555555555555556;
   GW(2) = .8888888888888889;
   GW(3) = .5555555555555556;
elseif (ng==4)
   GP(1) = -.86113631159405257524;
   GP(2) = -.33998104358485626481;
   GP(3) = .33998104358485626481;
   GP(4) = .86113631159405257524;
   GW(1) = .34785484513745385736;
   GW(2) = .65214515486254614264;
   GW(3) = .65214515486254614264;
   GW(4) = .34785484513745385736;
elseif (ng==5)
   GP(1) = -.90617984593866399282;
   GP(2) = -.53846931010568309105;
   GP(3) = 0;
   GP(4) = .53846931010568309105;
   GP(5) = .90617984593866399282;
   GW(1) = .23692688505618908749;
   GW(2) = .47862867049936646808;
   GW(3) = .56888888888888888888;
   GW(4) = .47862867049936646808;
   GW(5) = .23692688505618908749;
elseif (ng==6)
   GP(1) = -.9324695142031520;
   GP(2) = -.6612093864662645;
   GP(3) = -.2386191860831969;
   GP(4) = .2386191860831969;
   GP(5) = .6612093864662645;
   GP(6) = .9324695142031520;
   GW(1) = .1713244923791703;
   GW(2) = .3607615730481386;
   GW(3) = .4679139345726911;
   GW(4) = .4679139345726911;
   GW(5) = .3607615730481386;
   GW(6) = .1713244923791703;
elseif (ng==7)
   GP(1) = -.9491079123427585;
   GP(2) = -.7415311855993944;
   GP(3) = -.4058451513773972;
   GP(4) = 0;
   GP(5) = .4058451513773972;
   GP(6) = .7415311855993944;
   GP(7) = .9491079123427585;
   GW(1) = .1294849661688697;
   GW(2) = .2797053914892767;
   GW(3) = .3818300505051189;
   GW(4) = .4179591836734694;
   GW(5) = .3818300505051189;
   GW(6) = .2797053914892767;
   GW(7) = .1294849661688697;
elseif (ng==8)
   GP(1) = -.9602898564975362;
   GP(2) = -.7966664774136267;
   GP(3) = -.5255324099163290;
   GP(4) = -.1834346424956498;
   GP(5) = .1834346424956498;
   GP(6) = .5255324099163290;
   GP(7) = .7966664774136267;
   GP(8) = .9602898564975362;
   GW(1) = .1012285362903763;
   GW(2) = .2223810344533745;
   GW(3) = .3137066458778873;
   GW(4) = .3626837833783620;
   GW(5) = .3626837833783620;
   GW(6) = .3137066458778873;
   GW(7) = .2223810344533745;
   GW(8) = .1012285362903763;
elseif (ng==9)
   GP(1) = -.9681602395076261;
   GP(2) = -.8360311073266358;
   GP(3) = -.6133714327005904;
   GP(4) = -.3242534234038089;
   GP(5) = 0;
   GP(6) = .3242534234038089;
   GP(7) = .6133714327005904;
   GP(8) = .8360311073266358;
   GP(9) = .9681602395076261;
   GW(1) = .0812743883615744;
   GW(2) = .1806481606948574;
   GW(3) = .2606106964029354;
   GW(4) = .3123470770400029;
   GW(5) = .3302393550012598;
   GW(6) = .3123470770400028;
   GW(7) = .2606106964029355;
   GW(8) = .1806481606948574;
   GW(9) = .0812743883615744;
elseif (ng==10)
   GP(1) = -.973906528517172;
   GP(2) = -.865063366688985;
   GP(3) = -.679409568299024;
   GP(4) = -.433395394129247;
   GP(5) = -.148874338981631;
   GP(6) = .148874338981631;
   GP(7) = .433395394129247;
   GP(8) = .679409568299024;
   GP(9) = .865063366688985;
   GP(10) = .973906528517172;
   GW(1) = .066671344308688;
   GW(2) = .149451349150581;
   GW(3) = .219086362515982;
   GW(4) = .269266719309996;
   GW(5) = .295524224714753;
   GW(6) = .295524224714753;
   GW(7) = .269266719309996;
   GW(8) = .219086362515982;
   GW(9) = .149451349150581;
   GW(10) = .066671344308688;
elseif (ng==15)
   GP(1) = -.9879925180204854;
   GP(2) = -.9372733924007059;
   GP(3) = -.8482065834104272;
   GP(4) = -.7244177313601700;
   GP(5) = -.5709721726085388;
   GP(6) = -.3941513470775634;
   GP(7) = -.2011940939974345;
   GP(8) = 0;
   GP(9) = .2011940939974345;
   GP(10) = .3941513470775634;
   GP(11) = .5709721726085388;
   GP(12) = .7244177313601700;
   GP(13) = .8482065834104272;
   GP(14) = .9372733924007059;
   GP(15) = .9879925180204854;
   GW(1) = .03075324199611807;
   GW(2) = .07036604748811134;
   GW(3) = .1071592204671351;
   GW(4) = .1395706779261761;
   GW(5) = .1662692058169852;
   GW(6) = .1861610000155741;
   GW(7) = .1984314853271374;
   GW(8) = .2025782419255562;
   GW(9) = .1984314853271374;
   GW(10) = .1861610000155741;
   GW(11) = .1662692058169852;
   GW(12) = .1395706779261761;
   GW(13) = .1071592204671351;
   GW(14) = .07036604748811134;
   GW(15) = .03075324199611807;
end

return