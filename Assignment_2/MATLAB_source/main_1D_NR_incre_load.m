%main file to simulate 1D NBC loading in 1D linear elastic material
%incremental loading approach-based calculation is adopted here
clear; clc; close all;
%---------------(material properties)-------------------%
%---------------(material properties)-------------------%
E_0 = 100e6; % Young's modulus = 100 MPa
m = 40; % Material constant 'm' from problem statement
A = 1e-4; % Cross-sectional area

L = 1; % Length of bar = 1 m
nel = 10; % Number of FE elements
tolerance = 1e-8; % Relaxed slightly for non-linear convergence

outer_itrns = 101; % 100 steps of loading
inner_itern = 25; % NR inner iterations

%-----(Initial Tangent Stiffness Calculation)-----%
C_T = zeros(nel,1);
for i=1:nel
    % At initial state, strain is 0. 
    % The derivative of E*atan(m*e) evaluated at e=0 is E*m.
    C_T(i,1) = E_0 * m; 
end

%-----------------(Boundary condition)-----------------%
total_force = 10000; % Total applied force = 10 kN
traction_val = total_force / (outer_itrns - 1);  % Incremental load per step

del_u = 1;                          % DBC including (small disp at each outer iteration)
DBC_val_1st_node  = 0 ;                % All the case having zero disp at first node (fixed end)
NBC_val_last_node_C1 = traction_val ;  % NBC at last node
DBC_val_last_node_C2 = del_u ;         % DBC at last node

%--------------(2 Cases)---------------%
case_1 = 1;              % NBC
case_2 = 0;              % DBC
%--------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           ( real code start here )        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------(Nodal information  )------------------%
[nnp, lel, x_cord] = Nodal_information_fun(nel, L);
no_of_dof = nnp;
disp(' Nodal information calculated ');

%----------------( Mapping table)-------------------%
[element_nodes]= Mapping_table_fun(nel);
disp(' Mapping table genetrated');

%----------------(Stiffness matrix)-----------------%
%[K_global]= Stiffness_matrix_fun( nel, nnp, lel,  A, E,  x_cord,  element_nodes);
disp(' Global stiffness matrix genetrated ');

%======================================%
%  (NL_FEM starts here (NR_approach))  %
%======================================%
elemental_strain_history= zeros(outer_itrns, nel);  % strain preallocation
elemental_stress_history= zeros(outer_itrns, nel);  % stress preallocation
f_internal_history = zeros(outer_itrns, nnp);
U_global_history = zeros(outer_itrns, nnp) ; %storing deformation history
f_int_left = zeros(1, outer_itrns);%storing left node internal force
f_int_right = zeros(1, outer_itrns); %storing right node internal force
f_ext_right = zeros(1, outer_itrns);%storing applied force at right node
%---------------( disp array initialization)---------------%
% initialize all nodal disp ==0, DBC node will be modified as per DBC value
U_prev = zeros(nnp, 1); stress_prev = zeros(nel,1); strain_prev = zeros(nel,1);
f_int_prev = zeros(nnp,1);
tot_DBC_nodes = 0; index_nodes_DBC = zeros(1,2);

if (case_1==1 && case_2==0)                        % NBC
    BC=[DBC_val_1st_node, NBC_val_last_node_C1];
    tot_DBC_nodes = 1; index_nodes_DBC(1,1) = 1;
elseif (case_1==0 && case_2==1)                    % DBC
    BC=[DBC_val_1st_node, DBC_val_last_node_C2];   % DBC at left and right node
    tot_DBC_nodes = 2; index_nodes_DBC(1,1) = 1; index_nodes_DBC(1,2) = nnp;
end
F_ext_glbl = zeros(nnp, 1); F_ext_glbl(nnp,1) = traction_val;
f_int_global = zeros(nnp,1); f_int_cur = zeros(nnp,1); 
U_prev = zeros(nnp, 1);
u_cur_corr = zeros(nnp,1); f_ext_total = zeros(nnp,1);
%--------------------(outer disp iteration)--------------------%
for i = 2 : outer_itrns         %%%%(outer displacement loop)%%%%%
    f_ext_total = f_ext_total + F_ext_glbl;%total external traction applied so far
    disp('disp step = ');  disp(i);
    du_curr_total = zeros(nnp,1);%total incremental corrected disp in 
                                 % current outer iteration, which can also
                                 % be used to compute incremental strain
    %--------(jacobian)-------------------------%
    residual_cur = zeros(nnp,1);
    [Jacobian_global] = Comp_Jacobian(nel, nnp, lel, x_cord, element_nodes, C_T);
    [Jacobian_reduced, residual_reduced] = Jacobian_enforce_DBC_dbc_loading(no_of_dof, tot_DBC_nodes, index_nodes_DBC, ...
        Jacobian_global, F_ext_glbl, nnp, BC);
    %----------Jacobian computation + DBC----------
    du_cur_pred = (Jacobian_reduced)\residual_reduced;
    u_cur_pred = U_prev + du_cur_pred; %U_prev = total previous outer iteration displacement
    u_cur_corr = u_cur_pred; %total current displacement
    du_curr_total = du_cur_pred;

    [dstrn_cur_pred] = compute_strain(du_cur_pred, lel, nel);
    strn_cur_pred = dstrn_cur_pred + strain_prev;
    %--------update stress and internal force vector------------------------------
    [stress_cur_pred, C_T_all] = compute_VE_stress(strn_cur_pred, E_0, nel, strain_prev);
    [f_int_global] = compute_f_internal(stress_cur_pred, lel, nel, x_cord, element_nodes, nnp);
    %--------update stress and internal force vector------------------------------
    residual_cur = f_int_global - f_ext_total;

    residual_cur_free = zeros(nnp,1);
    residual_cur_free(2:(nnp),1)= residual_cur(2:(nnp),1); % corresponding to free DOFs
    
    %------inner iterations starting-----------
    if(ge(max(abs(residual_cur_free)), tolerance) == 1) %inner iteration
        du_corr_total = zeros(nnp,1); %total NR correction
        BC_tmp = zeros(1,2);
        for j=1:inner_itern
            [Jacobian_global] = Comp_Jacobian(nel, nnp, lel, x_cord, element_nodes, C_T);
            [Jacobian_reduced, residual_reduced] = Jacobian_enforce_DBC_dbc_loading(no_of_dof, tot_DBC_nodes, index_nodes_DBC, ...
                        Jacobian_global, residual_cur, nnp, BC_tmp);
            du_corr = - Jacobian_reduced \ residual_cur_free;
            du_corr_total = du_corr_total + du_corr;%incremental corrected displacement
            u_cur_corr = u_cur_pred + du_corr_total; %total corrected displacement
            du_curr_total = du_cur_pred + du_corr_total;%total incremental 
                                                        % displacement so far
            %du_curr_total: it can be used to compute incremental strain and 
            % add it to previous strain to get total current strain. 
            [dstrn_cur_pred] = compute_strain(du_curr_total, lel, nel);
            strn_cur_pred = dstrn_cur_pred + strain_prev;
            [stress_cur_pred, H_all] = compute_VE_stress(strn_cur_pred, E_0, nel, strain_prev);
            [f_int_global] = compute_f_internal(stress_cur_pred, lel, nel, x_cord, element_nodes, nnp);
            residual_cur = f_int_global - f_ext_total;
            residual_cur_free = zeros(nnp,1);
            residual_cur_free(2:(nnp),1)= residual_cur(2:(nnp),1); % corresponding to free DOFs
            if(le(max(abs(residual_cur_free)), tolerance) == 1)
                break;
            end
            if(j == inner_itern)
                disp('inner iterations count reached, there may be no convergence');
                disp('least residual =');
                max(abs(residual_cur_free))
       %         break;
            end
        end %for j=1:inner_iterns
    end %fi inner iterations
    %update previous converged iteration values in outer loop
    U_prev = u_cur_corr ; strain_prev = strn_cur_pred;
   f_ext_right(1,i) = f_ext_total(nnp,1);
    %---------update history variables---------------
    elemental_strain_history(i,:) = strn_cur_pred; elemental_stress_history(i,:) = stress_cur_pred;
    U_global_history(i,:) = u_cur_corr;
    f_internal_history(i,:) = f_int_global;
end %fi i=1:outer_iterns

% The following lines are added to plot the results after the simulation is complete.
% Figure 1: Stress-Strain curve for the last element
figure(1);
plot(elemental_strain_history(:, nel)*100, elemental_stress_history(:, nel)/10^6, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 12);
grid on;
xlabel('Strain (\epsilon) [%]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Stress (\sigma) [MPa]', 'FontSize', 12, 'FontWeight', 'bold');
title('Nonlinear Material Response (\sigma-\epsilon) in Last Element', 'FontSize', 14);
% Figure 2: Applied Force vs Displacement
figure(2);
plot(U_global_history(:, nnp), f_ext_right(1, :), 'r.-', 'LineWidth', 1.5, 'MarkerSize', 12);
grid on;
xlabel('Displacement at Right End Node (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Applied Force F (N)', 'FontSize', 12, 'FontWeight', 'bold');
title('Applied Force vs. Reaction Displacement', 'FontSize', 14);

% --- Save the plots for the LaTeX report ---
% Save Figure 1 (Stress-Strain) with 300 DPI resolution
exportgraphics(figure(1), 'stress_strain_plot.png', 'Resolution', 300);

% Save Figure 2 (Force-Displacement) with 300 DPI resolution
exportgraphics(figure(2), 'force_disp_plot.png', 'Resolution', 300);

disp('Plots successfully saved to the current directory!');