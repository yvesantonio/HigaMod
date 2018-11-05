function [AModified,bModified] = imposeBoundIGA(Al,bl,igaBoundCond,space,msh,coeff,geoData,spaceFunc,imb,kmb)

    %% NONHOMOGENEOUS DIRICHLET BOUNDARY CONDITIONS
    
    if(length(igaBoundCond.dirSides) == 2)
        
        disp('I am here!')
        
        % EXTRACT THE DIRICHLET INFORMATION
        
        [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, igaBoundCond.dir, igaBoundCond.dirSides);
        
        % COMPUTE THE LIFTING OPERATOR
        
        Lift  = @(x) (1 - x) * u_drchlt(1) + u_drchlt(2);
        dLift = @(x) -u_drchlt(1);
        
        % EVALUATE THE LIFTING OPERATOR IN THE HORIZONTAL QUADRATURE NODES
        
        evalLift = Lift(geoData.horNodes);
        evaldLift = dLift(geoData.horNodes);
        
        % COPUTE THE 1D BILINEAR FORM EVALUATED USING THE LIFTING OPERATOR
        
        numbKnots = msh.nel_dir;
        numbHorNodes = length(coeff.r00)/numbKnots;

        Local_00 = zeros (space.ndof, 1);
        Local_10 = zeros (space.ndof, 1);
        Local_01 = zeros (space.ndof, 1);
        Local_11 = zeros (space.ndof, 1);
        
        r00Lift = coeff.r00 .* evalLift;
        r10Lift = coeff.r10 .* evaldLift;
        r01Lift = coeff.r01 .* evalLift;
        r11Lift = coeff.r11 .* evaldLift;

        for iel = 1:numbKnots

            r00LiftLocal = r00Lift((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
            r10LiftLocal = r10Lift((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
            r01LiftLocal = r01Lift((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
            r11LiftLocal = r11Lift((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);

            Local_00 = Local_00 + op_f_v (spaceFunc{iel,1}, spaceFunc{iel,3}, r00LiftLocal);
            Local_10 = Local_10 + op_f_v (spaceFunc{iel,1}, spaceFunc{iel,3}, r10LiftLocal);
            Local_01 = Local_01 + op_f_v (spaceFunc{iel,1}, spaceFunc{iel,3}, r01LiftLocal);
            Local_11 = Local_11 + op_f_v (spaceFunc{iel,1}, spaceFunc{iel,3}, r11LiftLocal);

        end
        
        LiftRHS = Local_00 + Local_01 + Local_10 + Local_11;
        
        % IMPLEMENT MODIFICATIONS IN THE STIFFNESS MATRIX
        
        if (imb == kmb)
            AModified = Al;
            AModified(drchlt_dofs,:) = 0;
            AModified(drchlt_dofs,drchlt_dofs) = 1;
        else
            AModified = Al;
            AModified(drchlt_dofs,:) = 0;
        end
        
        % IMPLEMENT MODIFICATIONS IN THE FORCING TERM
        
        if (imb == kmb)
            bModified = bl + LiftRHS;
        else
            bModified = bl;
        end   
        
    elseif(length(igaBoundCond.dirSides) == 1)
        
        % EXTRACT THE DIRICHLET INFORMATION
        
        [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, igaBoundCond.dir, igaBoundCond.dirSides);
        
        % COMPUTE THE LIFTING OPERATOR
        
        Lift  = @(x) u_drchlt;
        dLift = @(x) 0;
        
        % EVALUATE THE LIFTING OPERATOR IN THE HORIZONTAL QUADRATURE NODES
        
        evalLift = Lift(geoData.horNodes);
        evaldLift = dLift(geoData.horNodes);
        
        % COPUTE THE 1D BILINEAR FORM EVALUATED USING THE LIFTING OPERATOR
        
        numbKnots = msh.nel_dir;
        numbHorNodes = length(coeff.r00)/numbKnots;

        Local_00 = zeros (space.ndof, 1);
        Local_10 = zeros (space.ndof, 1);
        Local_01 = zeros (space.ndof, 1);
        Local_11 = zeros (space.ndof, 1);
        
        r00Lift = coeff.r00 .* evalLift;
        r10Lift = coeff.r10 .* evaldLift;
        r01Lift = coeff.r01 .* evalLift;
        r11Lift = coeff.r11 .* evaldLift;

        for iel = 1:numbKnots

            r00LiftLocal = r00Lift((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
            r10LiftLocal = r10Lift((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
            r01LiftLocal = r01Lift((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);
            r11LiftLocal = r11Lift((iel - 1)*numbHorNodes + 1 : iel*numbHorNodes);

            Local_00 = Local_00 + op_f_v (spaceFunc{iel,1}, spaceFunc{iel,3}, r00LiftLocal);
            Local_10 = Local_10 + op_f_v (spaceFunc{iel,1}, spaceFunc{iel,3}, r10LiftLocal);
            Local_01 = Local_01 + op_f_v (spaceFunc{iel,2}, spaceFunc{iel,3}, r01LiftLocal);
            Local_11 = Local_11 + op_f_v (spaceFunc{iel,2}, spaceFunc{iel,3}, r11LiftLocal);

        end
        
        LiftRHS = Local_00 + Local_01 + Local_10 + Local_11;
        
        % IMPLEMENT MODIFICATIONS IN THE STIFFNESS MATRIX
        
        if (imb == kmb)
            AModified = Al;
            AModified(drchlt_dofs,:) = 0;
            AModified(drchlt_dofs,drchlt_dofs) = 1;
        else
            AModified = Al;
            AModified(drchlt_dofs,:) = 0;
        end
        
        % IMPLEMENT MODIFICATIONS IN THE FORCING TERM
        
        if (imb == kmb)
            bModified = bl + LiftRHS;
        else
            bModified = bl;
        end
        
    end
    
    %% NONHOMOGENEOUS ROBIN BOUNDARY CONDITIONS

end