function [evalDetJac] = detJac3D(x,y,z,Jac,type)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Change in input to take in consideration the specific
    %  parametrization of the 3D geometry
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xOld = x;
    yOld = y;
    zOld = z;
    
    if(strcmp(type,'Torus'))
        zNew = zOld;
        yNew = yOld;
        xNew = xOld;
    elseif(strcmp(type,'Cylinder'))
        zNew = xOld;
        yNew = yOld;
        xNew = zOld;
    end
    
    o = length(zNew);
    n = length(yNew);
    m = length(xNew);
    
    for ii = 1:n
        for jj = 1:m
            for kk = 1:o
                evalDetJac(ii,jj,kk)  = abs(det(Jac(xNew(jj),yNew(ii),zNew(kk))));
            end
        end
    end

end