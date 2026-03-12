%function to compute stress
function[stress_cur_pred, C_T_all] = compute_VE_stress(strn_cur_pred, E_0, nel, strain_prev)
strs_tmp = zeros(nel,1); C_T_all = zeros(nel,1);
for i=1:nel
    strs_tmp(i,1) = E_0 * strn_cur_pred(i,1); %compute appropriately
    C_T_all(i,1) = E_0; %compute appropriately
end
stress_cur_pred = strs_tmp;

