function   [Jacobian_global, Residual, D_cur]= Jacobian_residual_fun(U_prev, nel, nnp, lel,...
    x_cord, element_nodes, F_ext_glbl, History_variable, outer_itern)
Jacobian_local =zeros (2*nel,2) ;Jacobian_global = zeros(nnp, nnp); F_int_local=zeros( 2*nel ,1); F_int_glbl = zeros(nnp, 1);
%[elmn_sig_cur, elmn_strn, D_cur]= stress_update_fun(U_prev, lel, nel, W, E, m, D_tr);
[B] = Shape_function_fun(lel);     % Caling derivative of shape function

syms x
for i = 1:nel
    elmn_strn(i, 1)= B*U_prev([i,i+1],1);              % elemental strain
    
    elmn_stress(i,1)= (1-D_cur(i, 1))*E*elmn_strn(i, 1);% stress corresponding to newly updated damage-
    %-------------------------------
    dof_tmp= element_nodes(i, :) ;
    Jacobian_local([(2*i-1), 2*i],:) = int((B'* K_T(i, 1)* B), x, x_cord(i), x_cord(i+1)); % elememtal jacobian
    Jacobian_global(dof_tmp, dof_tmp) = Jacobian_global(dof_tmp, dof_tmp) + Jacobian_local([(2*i-1) , 2*i], :);
    F_int_local([(2*i-1), 2*i],:) =  int(B'*elmn_stress(i,1), x, x_cord(i), x_cord(i+1));  % elemental f_int
    F_int_glbl(dof_tmp, 1) = F_int_glbl( dof_tmp, 1) + F_int_local([(2*i-1), 2*i],:); % global f_int
end
F_ext_glbl = zeros(nnp, 1);      
Residual= F_int_glbl - F_ext_glbl;

%=========================================================================%


