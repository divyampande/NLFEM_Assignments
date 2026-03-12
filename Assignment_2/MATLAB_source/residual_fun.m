function [Residual, D_cur]=residual_fun(U_prev, lel, nel, W, E, x_cord, element_nodes, nnp, m, D_tr)
[B] = Shape_function_fun( lel);
F_int_local=zeros ( 2*nel ,1) ; F_int_glbl = zeros(nnp, 1);
[elmn_sig_cur, ~, D_cur]= stress_update_fun(U_prev, lel, nel, W, E, m, D_tr);
syms x                              % for integration-
for i = 1:nel
    F_int_local([(2*i-1), 2*i],:) =  int(B'*elmn_sig_cur(i,1), x, x_cord(i), x_cord(i+1));                 % elemental f_int
    F_int_glbl(element_nodes(i,:), 1) = F_int_glbl(element_nodes(i,:), 1) + F_int_local([(2*i-1), 2*i],:); % global f_int
end
F_ext_glbl = zeros(nnp, 1);        % will be updated for NBC case later.
Residual= F_int_glbl - F_ext_glbl;
