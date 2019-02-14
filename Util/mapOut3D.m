function [X,Y,Z] = mapOut3D(x,y,z,map,type)

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
    elseif(strcmp(type,'Slab_1'))
        zNew = xOld;
        yNew = yOld;
        xNew = zOld;
    elseif(strcmp(type,'Slab_2'))
        zNew = xOld;
        yNew = yOld;
        xNew = zOld;
    elseif(strcmp(type,'Slab_3'))
        zNew = xOld;
        yNew = yOld;
        xNew = zOld;
    elseif(strcmp(type,'Slab_4'))
        zNew = xOld;
        yNew = yOld;
        xNew = zOld;
    end
    
    l = length(zNew);
    n = length(yNew);
    m = length(xNew);
    X = zeros(m,n,l);
    Y = zeros(m,n,l);
    Z = zeros(m,n,l);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Change in output to take in consideration the specific
    %  parametrization of the 3D geometry
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ii = 1:m
        for jj = 1:n   
            for kk = 1:l
                mapped = map(xNew(ii),yNew(jj),zNew(kk));
                
                X(ii,jj,kk) = mapped(1);
                Y(ii,jj,kk) = mapped(2);
                Z(ii,jj,kk) = mapped(3);
                
            end
        end
    end
    
end