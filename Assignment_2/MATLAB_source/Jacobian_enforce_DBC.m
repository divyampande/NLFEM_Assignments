%function to enforce DBC on Jacobian matrix
function[Jacobian_reduced,residual_reduced] = Jacobian_enforce_DBC(no_of_dof, tot_DBC_nodes, index_nodes_DBC, ...
    Jacobian_global, residual_cur)
Jacobian_global_tmp = Jacobian_global; residual_tmp = residual_cur;
for i=1:tot_DBC_nodes
    Jacobian_global_tmp(index_nodes_DBC(1,i), 1:no_of_dof) = zeros(1,no_of_dof);
    Jacobian_global_tmp(1:no_of_dof,index_nodes_DBC(1,i)) = zeros(1, no_of_dof);
    Jacobian_global_tmp(index_nodes_DBC(i),index_nodes_DBC(i)) = 1;
    residual_tmp(index_nodes_DBC(i)) = 0;
end

Jacobian_reduced = Jacobian_global_tmp; residual_reduced = residual_tmp;