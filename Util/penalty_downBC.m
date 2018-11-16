function [r00Down, penalty_forceDown] = penalty_downBC(imb, kmb, mb_i_1, mb_yi_1, mb_k_1, mb_yk_1, Psi_y, LateralConditions)
    if strcmp(LateralConditions.bc_down, 'dir')
        penalty_forceDown =  LateralConditions.dato_down * mb_i_1;
        r00Down = mb_i_1 * mb_k_1;
    elseif strcmp(LateralConditions.bc_down, 'rob')
        penalty_forceDown  =  LateralConditions.dato_down * ( - LateralConditions.mu_down * Psi_y * mb_yi_1 + ...
                                                            LateralConditions.chi_down * mb_i_1 );
        r00Down = ( - LateralConditions.mu_down * Psi_y * mb_yi_1 + LateralConditions.chi_down * mb_i_1 ) .* ...
                  ( - LateralConditions.mu_down * Psi_y * mb_yk_1 + LateralConditions.chi_down * mb_k_1 );
    else
        disp('ATTENTION: unrecognized boundary condition!');
    end
end
