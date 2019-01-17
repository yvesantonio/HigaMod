function [evalJac,structPhi,evalDetJac] = jacOut3DHiMod(x,y,z,geoInfo,type)

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
    elseif(strcmp(type,'Slab'))
        zNew = xOld;
        yNew = zOld;
        xNew = yOld;
    end
    
    o = length(zNew);
    n = length(yNew);
    m = length(xNew);
    
    evalJac = cell(m,n,o);
    
    [Xref,Yref,Zref] = ndgrid(xNew,yNew,zNew);
    
    j11 = geoInfo.J11(Xref,Yref,Zref);
    j12 = geoInfo.J12(Xref,Yref,Zref);
    j13 = geoInfo.J13(Xref,Yref,Zref);
    j21 = geoInfo.J21(Xref,Yref,Zref);
    j22 = geoInfo.J22(Xref,Yref,Zref);
    j23 = geoInfo.J23(Xref,Yref,Zref);
    j31 = geoInfo.J31(Xref,Yref,Zref);
    j32 = geoInfo.J32(Xref,Yref,Zref);
    j33 = geoInfo.J33(Xref,Yref,Zref);
    
    evalDetJac = geoInfo.detJ(Xref,Yref,Zref);
    
    for ii = 1:m
        for jj = 1:n
            for kk = 1:o
                
                evalJac{ii,jj,kk} = [j11(ii,jj,kk) j12(ii,jj,kk) j13(ii,jj,kk);...
                                     j21(ii,jj,kk) j22(ii,jj,kk) j23(ii,jj,kk);...
                                     j31(ii,jj,kk) j32(ii,jj,kk) j33(ii,jj,kk)];
                                 
                structPhi.Phi1_dx(ii,jj,kk) = (j22(ii,jj,kk)*j33(ii,jj,kk) - ...
                                               j32(ii,jj,kk)*j23(ii,jj,kk))/ ...
                                               evalDetJac(ii,jj,kk);
                structPhi.Phi1_dy(ii,jj,kk) = (j13(ii,jj,kk)*j32(ii,jj,kk) - ...
                                               j33(ii,jj,kk)*j12(ii,jj,kk))/ ...
                                               evalDetJac(ii,jj,kk);
                structPhi.Phi1_dz(ii,jj,kk) = (j12(ii,jj,kk)*j23(ii,jj,kk) - ...
                                               j22(ii,jj,kk)*j13(ii,jj,kk))/ ...
                                               evalDetJac(ii,jj,kk);
                structPhi.Phi2_dx(ii,jj,kk) = (j23(ii,jj,kk)*j31(ii,jj,kk) - ...
                                               j33(ii,jj,kk)*j21(ii,jj,kk))/ ...
                                               evalDetJac(ii,jj,kk);
                structPhi.Phi2_dy(ii,jj,kk) = (j11(ii,jj,kk)*j33(ii,jj,kk) - ...
                                               j31(ii,jj,kk)*j13(ii,jj,kk))/ ...
                                               evalDetJac(ii,jj,kk);
                structPhi.Phi2_dz(ii,jj,kk) = (j13(ii,jj,kk)*j21(ii,jj,kk) - ...
                                               j23(ii,jj,kk)*j11(ii,jj,kk))/ ...
                                               evalDetJac(ii,jj,kk);
                structPhi.Phi3_dx(ii,jj,kk) = (j21(ii,jj,kk)*j32(ii,jj,kk) - ...
                                               j31(ii,jj,kk)*j22(ii,jj,kk))/ ...
                                               evalDetJac(ii,jj,kk);
                structPhi.Phi3_dy(ii,jj,kk) = (j12(ii,jj,kk)*j31(ii,jj,kk) - ...
                                               j32(ii,jj,kk)*j11(ii,jj,kk))/ ...
                                               evalDetJac(ii,jj,kk);
                structPhi.Phi3_dz(ii,jj,kk) = (j11(ii,jj,kk)*j22(ii,jj,kk) - ...
                                               j21(ii,jj,kk)*j12(ii,jj,kk))/ ...
                                               evalDetJac(ii,jj,kk);
                
            end
        end
    end

end