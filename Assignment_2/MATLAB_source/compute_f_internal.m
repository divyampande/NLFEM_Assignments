%function to compute internal force vector
function[f_int_global] = compute_f_internal(stress_cur_pred, lel, nel, x_cord, element_nodes, nnp)
[B] = Shape_function_fun( lel);
F_int_local=zeros ( 2*nel ,1) ; F_int_glbl = zeros(nnp, 1);
syms x                              % for integration-
for i = 1:nel
    F_int_local([(2*i-1), 2*i],:) =  int(B'*stress_cur_pred(i,1), x, x_cord(i), x_cord(i+1));                 % elemental f_int
    F_int_glbl(element_nodes(i,:), 1) = F_int_glbl(element_nodes(i,:), 1) + F_int_local([(2*i-1), 2*i],:); % global f_int
end

f_int_global = F_int_glbl;