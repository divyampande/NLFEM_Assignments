%function to enforce DBC on Jacobian matrix when DBC condition at the end
function[Jacobian_reduced,residual_reduced] = Jacobian_enforce_DBC_dbc_loading(no_of_dof, tot_DBC_nodes, index_nodes_DBC, ...
    Jacobian_global, residual_cur, nnp, BC)
Jacobian_global_tmp = Jacobian_global; residual_tmp = residual_cur;
tmp2 = zeros(no_of_dof,1);
for i=1:no_of_dof
    tmp = Jacobian_global_tmp(i,index_nodes_DBC(1:tot_DBC_nodes)).*BC(1:tot_DBC_nodes);
    tmp = sum(tmp);
    tmp2(i,1) = tmp;
end

residual_tmp = residual_tmp - tmp2;
residual_tmp(index_nodes_DBC(1:tot_DBC_nodes),1) = BC(1:tot_DBC_nodes);

for i=1:tot_DBC_nodes
    Jacobian_global_tmp(index_nodes_DBC(1,i), 1:no_of_dof) = zeros(1,no_of_dof);
    Jacobian_global_tmp(1:no_of_dof,index_nodes_DBC(1,i)) = zeros(1, no_of_dof);
    Jacobian_global_tmp(index_nodes_DBC(i),index_nodes_DBC(i)) = 1;
end

Jacobian_reduced = Jacobian_global_tmp; residual_reduced = residual_tmp;