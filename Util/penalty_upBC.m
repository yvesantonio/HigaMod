function [r00Up, penalty_forceUp] = penalty_upBC(imb, kmb, mb_i_end, mb_yi_end, mb_k_end, mb_yk_end, Psi_y, LateralConditions)
    if strcmp(LateralConditions.bc_up, 'dir')
        penalty_forceUp  =  LateralConditions.dato_up * mb_i_end;
        r00Up = mb_i_end * mb_k_end;
    elseif strcmp(LateralConditions.bc_up, 'rob')
        penalty_forceUp  =  LateralConditions.dato_up * ( LateralConditions.mu_up * Psi_y * mb_yi_end + ...
                                                          LateralConditions.chi_up * mb_i_end );
        r00Up = ( LateralConditions.mu_up * Psi_y * mb_yi_end + LateralConditions.chi_up * mb_i_end ) .* ...
                ( LateralConditions.mu_up * Psi_y * mb_yk_end + LateralConditions.chi_up * mb_k_end );
    else
        disp('ATTENTION: unrecognized boundary condition!');
    end
end

