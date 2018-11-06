function [structPhi] = invJac3D(evalJac)

    [n,m,o] = size(evalJac);
    
    for ii = 1:n
        for jj = 1:m
            for kk = 1:o
                evalJavInv = evalJac{ii,jj,kk}\eye(size(evalJac{ii,jj,kk}));
                
                structPhi.Phi1_dx(ii,jj) = evalJavInv(1,1);
                structPhi.Phi1_dy(ii,jj) = evalJavInv(1,2);
                structPhi.Phi1_dz(ii,jj) = evalJavInv(1,3);
                structPhi.Phi2_dx(ii,jj) = evalJavInv(2,1);
                structPhi.Phi2_dy(ii,jj) = evalJavInv(2,2);
                structPhi.Phi2_dz(ii,jj) = evalJavInv(2,3);
                structPhi.Phi3_dx(ii,jj) = evalJavInv(3,1);
                structPhi.Phi3_dy(ii,jj) = evalJavInv(3,2);
                structPhi.Phi3_dz(ii,jj) = evalJavInv(3,3);
            end
        end
    end

end