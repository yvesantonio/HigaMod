function [AlIn,blIn] = penalty_inflowBC(imb, kmb, augVerWeights, mb_i, mb_k, detJac, ComputedInOut, BoundCond, space, spaceFunc, msh)
    AlIn = spalloc (space.ndof, space.ndof, 3*space.ndof);
    blIn = zeros (space.ndof, 1);
    
    iel = 1;
    dtheta = zeros(space.ndof,1); 
    theta  = zeros(space.ndof,1);    
    dtheta(1:3) = dtheta_IO(spaceFunc{iel,2}, spaceFunc{iel,3}, 0);
    theta(1:3)  = theta_IO(spaceFunc{iel,1}, spaceFunc{iel,3}, 0);
    
    % mettiamo a zero perché sì
%     dtheta(3) = 0; 
%     theta(1:3)  = [1 0 0];
    
    if BoundCond.dirSides(1) == 1
       AlIn(1,1) = sum( detJac .* mb_i .* mb_k .* augVerWeights );
       blIn(1)   = sum( detJac .* ComputedInOut.dirIn_c .* mb_i .* augVerWeights );
    elseif BoundCond.robSides(1) == 1
        for j=1:length(blIn)
            for l=1:length(blIn)                
                AlIn(j,l) = sum( ( -dtheta(l) * ComputedInOut.muIn_c + ComputedInOut.chiIn_c .* theta(l) ) .* ...
                                 ( -dtheta(j) * ComputedInOut.muIn_c + ComputedInOut.chiIn_c .* theta(j) ) .* ... 
                                  mb_i .* mb_k .* detJac .* augVerWeights );
            end
            blIn(j) = sum( ( -dtheta(j) * ComputedInOut.muIn_c + ComputedInOut.chiIn_c .* theta(j) ) .* ...
                           ComputedInOut.robIn_c .* mb_i .* detJac .* augVerWeights );
        end
    end
end

