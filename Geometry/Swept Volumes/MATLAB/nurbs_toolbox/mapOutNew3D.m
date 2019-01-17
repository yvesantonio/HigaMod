function [varIn,Xout,Yout,Zout,selectOUT] = mapOutNew3D(x,y,z,map,out,type)

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
    
    l = length(zNew);
    n = length(yNew);
    m = length(xNew);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Change in output to take in consideration the specific
    %  parametrization of the 3D geometry
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(strcmp(type,'Torus'))
        if(out == 1) 
            outNew = 1;
        elseif(out == 2) 
            outNew = 3;
        elseif(out == 3) 
            outNew = 2;
        end
    elseif(strcmp(type,'Cylinder'))
        if(out == 1) 
            outNew = 1;
        elseif(out == 2) 
            outNew = 2;
        elseif(out == 3) 
            outNew = 3;
        end
    end
    
    [XX,YY,ZZ] = meshgrid(xNew,yNew,zNew);
    
    X = reshape(XX,[],1);
    Y = reshape(YY,[],1);
    Z = reshape(ZZ,[],1);
    
    varIn = [X,Y,Z];
    
    mapped = map(X',Y',Z');
    
    Xout = mapped(1,:);
    Yout = mapped(1,:);
    Zout = mapped(1,:);
    
    XOUT = reshape(Xout,n,m,l);
    YOUT = reshape(Yout,n,m,l);
    ZOUT = reshape(Zout,n,m,l);
    
    if (out == 1) 
        selectOUT = XOUT;
    elseif (out == 2)
        selectOUT = YOUT;
    else
        selectOUT = ZOUT;
    end
    
end