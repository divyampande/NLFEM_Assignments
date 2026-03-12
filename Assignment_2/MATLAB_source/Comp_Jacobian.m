%function to compute Jacobian matrix for VE case
function[Jacobian_global]= Comp_Jacobian(nel, nnp, lel,...
    x_cord, element_nodes, C_T)
Jacobian_local =zeros (2*nel,2) ;Jacobian_global = zeros(nnp, nnp); 
[B] = Shape_function_fun(lel);     % Caling derivative of shape function

syms x
for i = 1:nel
    dof_tmp= element_nodes(i, :) ;
    Jacobian_local([(2*i-1), 2*i],:) = int((B'* C_T(i,1) * B), x, x_cord(i), x_cord(i+1)); % elememtal jacobian
    Jacobian_global(dof_tmp, dof_tmp) = Jacobian_global(dof_tmp, dof_tmp) + Jacobian_local([(2*i-1) , 2*i], :);
end
%=========================================================================%


