%function to compute stress
function [stress_cur_pred, C_T_all] = compute_VE_stress(strn_cur_pred, E_0, nel, ~)
    strs_tmp = zeros(nel, 1); 
    C_T_all = zeros(nel, 1);
    
    m = 40; % Material constant from problem statement
    
    for i = 1:nel
        epsilon = strn_cur_pred(i, 1); 
        
        % Non-linear stress equation
        strs_tmp(i, 1) = E_0 * atan(m * epsilon); 
        
        % Tangent Modulus (derivative of stress wrt strain)
        C_T_all(i, 1) = (E_0 * m) / (1 + (m * epsilon)^2); 
    end
    
    stress_cur_pred = strs_tmp;
end