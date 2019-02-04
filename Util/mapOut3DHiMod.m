function [X,Y,Z] = mapOut3DHiMod(x,y,z,geoInfo,type)

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
        yNew = zOld;
        xNew = yOld;
    elseif(strcmp(type,'Slab_2'))
        zNew = xOld;
        yNew = zOld;
        xNew = yOld;
    elseif(strcmp(type,'Slab_3'))
        zNew = xOld;
        yNew = zOld;
        xNew = yOld;
    end
    
    l = length(zNew);
    n = length(yNew);
    m = length(xNew);
    X = zeros(n,m,l);
    Y = zeros(n,m,l);
    Z = zeros(n,m,l);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Change in output to take in consideration the specific
    %  parametrization of the 3D geometry
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Xref,Yref,Zref] = ndgrid(xNew,yNew,zNew);
            
    X = geoInfo.Phi1(Xref,Yref,Zref);
    Y = geoInfo.Phi2(Xref,Yref,Zref);
    Z = geoInfo.Phi3(Xref,Yref,Zref);
    
%     for ii = 1:n
%         for jj = 1:m   
%             for kk = 1:l
%                 mapped = map(xNew(jj),yNew(ii),zNew(kk));
%                 
%                 X(ii,jj,kk) = mapped(1);
%                 Y(ii,jj,kk) = mapped(2);
%                 Z(ii,jj,kk) = mapped(3);
%                 
%             end
%         end
%     end
    
end