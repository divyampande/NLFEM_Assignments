function [ Internal_nodes ]= BC_imposition_fun(BC, nnp, case_1, case_2 )

%---------(Coloumn vector of known dof)---------------------------%
U_known_dof  = zeros(nnp,1);  % Preallocation for known dof vector
FB_global= zeros(nnp,1);      % Preallocation for global FB vector % unknown

if case_1==1 && case_2==0  % NBC

   % U_known_dof(1 ,1) = BC(1);                  % dirichlet BC_1
    Internal_nodes =  2 : nnp ;                 % Removing row and coloumn corresponding to DBC

%     Tmp_vec = K_global * U_known_dof ;          % To capture the effect of non zero dbc
% 
%     FB_global(nnp,1) =  A * BC(2) ;
%     F_global = FI_global + FB_global;
%     F_global= F_global - Tmp_vec;               % Modified global forced vector for row and coloumn elimination-

elseif case_1==0 && case_2==1  % DBC
   
    U_known_dof(1 ,1) = BC(1);                  % dirichlet BC_1
    U_known_dof(nnp ,1) = BC(2);                % dirichlet BC_2
    
    Internal_nodes =  2:(nnp-1);                % Removing row and coloumn corresponding to DBC

%     Tmp_vec = K_global * U_known_dof;           % To capture the effect of non zero dbc
% 
%     F_global = FI_global + FB_global;
%     F_global = F_global - Tmp_vec;              % Modified global forced vector for row and coloumn elimination-
end


end