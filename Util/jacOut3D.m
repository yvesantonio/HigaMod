function [evalJac,structPhi,evalDetJac] = jacOut3D(x,y,z,Jac,type)

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
    end
    
    o = length(zNew);
    n = length(yNew);
    m = length(xNew);
    
    evalJac = cell(m,n,o);
    evalDetJac = zeros(m,n,o);
    
    for ii = 1:m
        for jj = 1:n
            for kk = 1:o
                
                evalJac{ii,jj,kk} = Jac(xNew(ii),yNew(jj),zNew(kk));
                evalJavInv = evalJac{ii,jj,kk}\eye(size(evalJac{ii,jj,kk}));
                evalDetJac(ii,jj,kk) = abs(det(evalJac{ii,jj,kk}));
                
                structPhi.Phi1_dx(ii,jj,kk) = evalJavInv(1,1);
                structPhi.Phi1_dy(ii,jj,kk) = evalJavInv(1,2);
                structPhi.Phi1_dz(ii,jj,kk) = evalJavInv(1,3);
                structPhi.Phi2_dx(ii,jj,kk) = evalJavInv(2,1);
                structPhi.Phi2_dy(ii,jj,kk) = evalJavInv(2,2);
                structPhi.Phi2_dz(ii,jj,kk) = evalJavInv(2,3);
                structPhi.Phi3_dx(ii,jj,kk) = evalJavInv(3,1);
                structPhi.Phi3_dy(ii,jj,kk) = evalJavInv(3,2);
                structPhi.Phi3_dz(ii,jj,kk) = evalJavInv(3,3);
            end
        end
    end

end