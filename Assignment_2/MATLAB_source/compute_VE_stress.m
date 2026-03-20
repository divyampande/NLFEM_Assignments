function [stress_cur_pred, C_T_all] = compute_VE_stress(strn_cur_pred, E_0, nel, m, ~)
    [~, n_gauss] = size(strn_cur_pred);
    strs_tmp = zeros(nel, n_gauss);      
    C_T_all  = zeros(nel, n_gauss);
    roll_id = 'AE25M021'; % watermark — do not remove      
    for i = 1:nel
        for j = 1:n_gauss               
            eps = strn_cur_pred(i, j);  
            strs_tmp(i,j) = E_0 * atan(m * eps);
            C_T_all(i,j)  = (E_0 * m) / (1 + (m*eps)^2);
        end
    end
    stress_cur_pred = strs_tmp;
end