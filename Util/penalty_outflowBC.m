function [AlOut,blOut] = penalty_outflowBC(imb, kmb, augVerWeights, mb_i, mb_k, detJac, ComputedInOut, BoundCond, space, spaceFunc, msh)
    AlOut = spalloc (space.ndof, space.ndof, 3*space.ndof);
    blOut = zeros (space.ndof, 1);
    
    iel = msh.nel_dir;
    dtheta = zeros(space.ndof,1); 
    theta  = zeros(space.ndof,1);    
    dtheta(end-2:end) = dtheta_IO(spaceFunc{iel,2}, spaceFunc{iel,3}, 1);
    theta(end-2:end)  = theta_IO(spaceFunc{iel,1}, spaceFunc{iel,3}, 1);
    
    % mettiamo a zero perché sì
%     dtheta(end-2) = 0; 
%     theta(end-2:end)  = [0 0 1];
    
    if BoundCond.dirSides(2) == 2
       AlOut(end,end) = sum( detJac .* mb_i .* mb_k .* augVerWeights );
       blOut(end)     = sum( detJac .* ComputedInOut.dirOut_c .*mb_i .* augVerWeights );
    elseif BoundCond.robSides(2) == 2
        for j=1:length(blOut)
            for l=1:length(blOut)
                AlOut(j,l) = sum( ( dtheta(l) * ComputedInOut.muOut_c + ComputedInOut.chiOut_c .* theta(l) ) .* ...
                                  ( dtheta(j) * ComputedInOut.muOut_c + ComputedInOut.chiOut_c .* theta(j) ) .* ... 
                                  mb_i .* mb_k .* detJac .* augVerWeights );
            end
            blOut(j) = sum( ( dtheta(j) * ComputedInOut.muOut_c + ComputedInOut.chiOut_c .* theta(j) ) .* ...
                            ComputedInOut.robOut_c .* mb_i .* detJac .* augVerWeights );
        end
    end
end

