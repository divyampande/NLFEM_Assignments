function [U]=Solution_fun(BC, nnp, K_global, Internal_nodes, F_global, case_1, case_2 )

U= zeros(nnp,1);     % Preallocation of displacement vector

%----------(Reduced stiffness matrix and forced vector)------------%

K_red = K_global(Internal_nodes, Internal_nodes);
F_red = F_global(Internal_nodes,1);

%----------(Displacement vector solution)-----------------------------%
if case_1==1 && case_2==0      % FTET
    U(1,1) = BC(1);            % dirichlet_BC_1
   
elseif case_1==0 && case_2==1  % RTGT

    U(1,1)  = BC(1);          % dirichlet_BC_1
    U(nnp,1)= BC(2);          % dirichlet_BC_2

end

U_red = inv(K_red) * F_red ;      % Calculation of reduced displacement vector      
U(Internal_nodes, 1)= U_red;

end 