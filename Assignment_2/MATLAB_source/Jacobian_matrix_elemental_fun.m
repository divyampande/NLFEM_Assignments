function   [Jacobian_global]= jacobian_matrix_elemental_fun(U_prev, nel, nnp, lel,...
    x_cord, element_nodes, D_prev, E, m, W)
Jacobian_local =zeros (2*nel,2) ;Jacobian_global = zeros(nnp, nnp); F_int_local=zeros( 2*nel ,1); F_int_glbl = zeros(nnp, 1);
%[elmn_sig_cur, elmn_strn, D_cur]= stress_update_fun(U_prev, lel, nel, W, E, m, D_tr);
[B] = Shape_function_fun( lel);     % Caling derivative of shape function

syms x
for i = 1:nel
    elmn_strn(i, 1)= B*U_prev([i,i+1],1);              % elemental strain
    elmn_strn_enrgy(i, 1)= (1/2)*E*elmn_strn(i)^2;     % elemental strain energy
    f_tr(i, 1)=(1-D_prev(i, 1))^(1/ m)*elmn_strn_enrgy(i, 1)- W;
    %------------------------------
    if (le(f_tr(i, 1),0))    % if point is under surface-
        D_cur= D_prev; 
        K_T(i, 1)= E;
    else
        D_cur(i, 1)= 1- (W/elmn_strn_enrgy(i, 1))^(1/m);
        K_T(i, 1)= (1-D_cur(i, 1))*(1-(2/m))*E;         % point is on or out of energy surface
    end
    
    %-------------------------------
    dof_tmp= element_nodes(i, :) ;
    Jacobian_local([(2*i-1), 2*i],:) = int((B'* K_T(i, 1)* B), x, x_cord(i), x_cord(i+1)); % elememtal jacobian
    Jacobian_global(dof_tmp, dof_tmp) = Jacobian_global(dof_tmp, dof_tmp) + Jacobian_local([(2*i-1) , 2*i], :);
   
end

%=========================================================================%


