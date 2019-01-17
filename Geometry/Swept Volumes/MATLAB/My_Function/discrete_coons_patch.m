function Coons_patch = discrete_coons_patch(c1, c2, c3, c4, nbf_u, nbf_v)
% Function to compute the control points of the Coons patch of a given contour
% defined by four curves, c1, c2, c3, and c4
%
% (c) Nicolas Douillet, 2017
%
%%% Syntax
% Coons_patch = discrete_coons_patch(c1, c2, c3, c4)
% Coons_patch = discrete_coons_patch(c1, c2, c3, c4, nbf_u, nbf_v)
%
%%% Description
% Coons_patch = discrete_coons_patch(c1, c2, c3, c4) compute the Coons
% patch of the quadrangle (c1, c2, c3, c4) with nbf_u = nbf_v = 5
% (functions) by default
%
% Coons_patch = discrete_coons_patch(c1, c2, c3, c4, nbf_u, nbf_v) uses
% nbf_u functions in u direction, and nbf_v functions in v direction.
%
%%% Curves orientation
%
%              ^ v
%              |
%              |                c3
%        (1,0) o- - - - - - - - > - - - - - - - -o (1,1)
%              |                                 |
%              |                                 |
%              |                                 |
%           c2 ^              S                  ^ c4
%              |                                 |
%              |                                 |
%              |                                 |
%              o- - - - - - - - > - - - - - - - -o- - - > u
%          (0,0)                c1               (0,1)
%
% Connection constraints in X, Y, and Z :
%
% c1(1) = c3(1)
% c1(end) = c4(1)
% c2(1) = c3(end)
% c2(end) = c4(end)
%
%% Input arguments %%
%
% - c1 : curve #1, alongside u parameter (see accepted formats below)
% - c2 : curve #2, alongside v parameter (see accepted formats below)
% - c3 : curve #3, alongside u parameter, on the opposite side of c1
% - c4 : curve #4, alongside v parameter, on the opposite side of c2
% - nbf_u : number of functions / knots in u direction (integer >= 3)
% - nbf_v : number of functions / knots in v direction (integer >= 3)
%
% (c1, c2, c3, c4) / cn accepted formats
%         |       |       |
% - cn = [X] ... [Y] ... [Z] => size(cn) = [nb_sample 1 3] (RECOMMENDED FORMAT)
%         |       |       |
%
% or
%
%        | | |
% - cn = X Y Z => size(cn) = [nb_sample 3] or
%        | | |
%
%        -- X --
% - cn = -- Y -- => size(cn) = [3 nb_sample]
%        -- Z --
%
%%% Sampling advice for C1, c2, c3, and c4
% Take a sufficient (>~ 10) samples / number of elements for c1, c2, c3,
% and c4 vector regarding the numbers nbf_u and nbf_v of functions choosen
% in u and v direction. Mind also the Shannon Nyquist sampling condition.
%
%% Output argument %%
%
% - Coons_patch : a matrix / array of size [nbf_u, nbf_v, 3] giving
%   the NURBS control points coordinates
% 
%   such that     X = Coons_patch(:,:,1)
%                 Y = Coons_patch(:,:,2)
%                 Z = Coons_patch(:,:,3)
%
%% Example %%
%
%%% 3D parabole example and plot the result
% nb_sample = 64;
% x_min = 0;
% x_max = 9;
% y_min = 0;
% y_max = 9
% 
% c1 = cat(3, cat(3, linspace(x_min,x_max,nb_sample)', y_min*ones(nb_sample,1)), (linspace(x_min,x_max,nb_sample)').^2);
% c3 = cat(3, cat(3, c1(:,1,1), y_max*ones(nb_sample,1), (y_max*ones(nb_sample,1)).^2 + c1(:,1,3)));
% c2 = cat(3, cat(3, c1(:,1,2), c1(:,1,1)), c1(:,1,3));
% c4 = cat(3, cat(3, c3(:,1,2), c2(:,1,2)), c3(:,1,3));
% 
% % Surface
% S = zeros(nb_sample, nb_sample, 3);
% [S(:,:,1) S(:,:,2)] = meshgrid(linspace(x_min, x_max, nb_sample), linspace(x_min, x_max, nb_sample));
% S(:,:,3) = S(:,:,1).^2 + S(:,:,2).^2;
% 
% nbf_u = 8;
% nbf_v = 8;
% 
% P = discrete_coons_patch(c1, c2, c3, c4, nbf_u, nbf_v);
%
% figure;
% % Contour : C1, c2, c3, and c4 curves 
% line(c1(:,1,1),c1(:,1,2),c1(:,1,3), 'Color', [1 0 0], 'Linewidth',2), hold on;
% line(c2(:,1,1),c2(:,1,2),c2(:,1,3), 'Color', [1 0 0], 'Linewidth',2), hold on;
% line(c3(:,1,1),c3(:,1,2),c3(:,1,3), 'Color', [1 0 0], 'Linewidth',2), hold on;
% line(c4(:,1,1),c4(:,1,2),c4(:,1,3), 'Color', [1 0 0], 'Linewidth',2), hold on;
% % Interpolated -assumed- surface 
% surf(S(:,:,1), S(:,:,2), S(:,:,3)), shading interp, hold on;
% % Coons patch
% plot3(P(1:end,1:end,1),P(1:end,1:end,2),P(1:end,1:end,3),'ko','Linewidth',2), hold on;
% for i = 2:nbf_u-1
%     plot3(P(i,1:end,1),P(i,1:end,2),P(i,1:end,3),'Color','k','Linewidth',2), hold on;
% end
% for j = 2:nbf_v-1
%     plot3(P(1:end,j,1),P(1:end,j,2),P(1:end,j,3),'Color','k','Linewidth',2), hold on;
% end
% axis square;
% view(3);
%
%% Input Wrangling %%
%
if nargin < 6
   nbf_v = 5; 
end

if nargin < 5
   nbf_u = 5;
end

assert(nargin >= 4, 'Not enough input arguments');
assert(nbf_u >= 3 && isinteger(uint8(nbf_u)), 'Number of functions in u direction must be an integer greater or equal to 3.');
assert(nbf_v >= 3 && isinteger(uint8(nbf_v)), 'Number of functions in v direction must be an integer greater or equal to 3.');

% Curves format check and reshape
assert(numel(c1) >= 12 && ...
       numel(c2) >= 12 && ...
       numel(c3) >= 12 && ...
       numel(c4) >= 12, ...
       'Curves c1, Cé, c3, and c4 must contain at least twelve data -four 3D points- each.');

s1 = size(c1);
s2 = size(c2);
s3 = size(c3);
s4 = size(c4);

assert(s1(1) == 3 || s1(2) == 3 || s1(3) == 3, 'c1 : bad data format. Type ''help discrete_coons_patch'' to see the list of accepted formats.');

if(s1(3) == 1)
    if(s1(1) == 3 && s1(2) > 3)
        c1 = c1';
    end
    c1 = reshape(c1, [s1(1) 1 3]);
end

if(s2(3) == 1)
    if(s2(1) == 3 && s2(2) > 3)
        c2 = c2';
    end
    c2 = reshape(c2, [s2(1) 1 3]);
end

if(s3(3) == 1)
    if(s3(1) == 3 && s3(2) > 3)
        c3 = c3';
    end
    c3 = reshape(c3, [s3(1) 1 3]);
end

if(s4(3) == 1)
    if(s4(1) == 3 && s4(2) > 3)
        c4 = c4';
    end
    c4 = reshape(c4, [s4(1) 1 3]);
end
%
%% Contour initialization %%
%
Coons_patch = zeros(nbf_u,nbf_v,3);

c1_spl_vect = zeros(nbf_u,1);
c1_spl_vect(1,1) = 1;
c1_spl_vect(end,1) = s1(1);
c1_spl_vect(2:end-1,1) = round((s1(1)) * round(100*(1:nbf_u-2)/(nbf_u-1))/100);

c3_spl_vect = zeros(nbf_u,1);
c3_spl_vect(1,1) = 1;
c3_spl_vect(end,1) = s3(1);
c3_spl_vect(2:end-1,1) = round((s3(1)) * round(100*(1:nbf_u-2)/(nbf_u-1))/100);

c2_spl_vect = zeros(nbf_v,1);
c2_spl_vect(1,1) = 1;
c2_spl_vect(end,1) = s2(1);
c2_spl_vect(2:end-1,1) = round((s2(1)) * round(100*(1:nbf_v-2)/(nbf_v-1))/100);

c4_spl_vect = zeros(nbf_v,1);
c4_spl_vect(1,1) = 1;
c4_spl_vect(end,1) = s4(1);
c4_spl_vect(2:end-1,1) = round((s4(1)) * round(100*(1:nbf_v-2)/(nbf_v-1))/100);

Coons_patch(:,1,:)   = c1(c1_spl_vect,1,:);
Coons_patch(:,end,:) = c3(c3_spl_vect,1,:);
Coons_patch(1,:,:)   = c2(c2_spl_vect,1,:);
Coons_patch(end,:,:) = c4(c4_spl_vect,1,:);
%
%% Patch values computation %%
%
for n = 1:nbf_u
    for p = 1:nbf_v

        u = (n-1)/(nbf_u-1);
        v = (p-1)/(nbf_v-1);

        Coons_patch(n,p,1) = (1-u)*Coons_patch(1,p,1)+...
                             (1-v)*Coons_patch(n,1,1)+...
                                 u*Coons_patch(end,p,1) +...
                                 v*Coons_patch(n,end,1) +...
                                 ...
                       (u-1)*(1-v)*Coons_patch(1,1,1)-...
                               u*v*Coons_patch(end,end,1)+...
                           u*(v-1)*Coons_patch(end,1,1)+...
                           v*(u-1)*Coons_patch(1,end,1);

        Coons_patch(n,p,2) = (1-u)*Coons_patch(1,p,2)+...
                             (1-v)*Coons_patch(n,1,2)+...
                                 u*Coons_patch(end,p,2) +...
                                 v*Coons_patch(n,end,2) +...
                                 ...
                       (u-1)*(1-v)*Coons_patch(1,1,2)-...
                               u*v*Coons_patch(end,end,2)+...
                           u*(v-1)*Coons_patch(end,1,2)+...
                           v*(u-1)*Coons_patch(1,end,2);

        Coons_patch(n,p,3) = (1-u)*Coons_patch(1,p,3)+...
                             (1-v)*Coons_patch(n,1,3)+...
                                 u*Coons_patch(end,p,3) +...
                                 v*Coons_patch(n,end,3) +...
                                 ...
                       (u-1)*(1-v)*Coons_patch(1,1,3)-...
                               u*v*Coons_patch(end,end,3)+...
                           u*(v-1)*Coons_patch(end,1,3)+...
                           v*(u-1)*Coons_patch(1,end,3);

    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%discrete_coons_patch

