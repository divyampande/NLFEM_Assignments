%function to compute total current strain
function[strn_cur_pred] = compute_strain(u_cur_pred, lel, nel)
[B] = Shape_function_fun(lel);
elmn_strn= zeros(nel,1);
for i = 1:nel
    elmn_strn(i, 1)= B*u_cur_pred([i,i+1],1);                       % elemental strain
end

strn_cur_pred = elmn_strn;