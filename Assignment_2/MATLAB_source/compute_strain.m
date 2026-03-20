function [strn_cur_pred] = compute_strain(u_cur_pred, lel, nel, n_gauss)
    [B] = Shape_function_fun(lel);
    elmn_strn = zeros(nel, n_gauss);
    for i = 1:nel
        u_e = u_cur_pred([i, i+1], 1);
        for j = 1:n_gauss              
            elmn_strn(i, j) = B * u_e; % B constant, same value at all GPs
        end
    end
    strn_cur_pred = elmn_strn;
end