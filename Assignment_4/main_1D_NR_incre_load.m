%main file to simulate 1D bar with elasto-plastic / elasto-viscoplastic response
%incremental loading approach-based calculation is adopted here
% clear; clc; close all;
%---------------(material properties)-------------------%
E        = 210e9;       % Young's modulus = 210 GPa (steel)
sigma_y  = 250e6;       % Initial yield stress = 250 MPa
H_iso    = 10e9;        % Bilinear isotropic hardening modulus
K_hard   = 500e6;       % Exponential hardening: saturated increment (Pa)
b_hard   = 15;          % Exponential hardening: saturation rate
eta_vp   = 1e4;         % Viscosity (Pa.s) — EVP only
xi_vp    = 0.5;           % Viscoplastic rate exponent — EVP only
A        = 1e-4;        % Cross-sectional area (m^2)

L          = 1;         % Length of bar = 1 m
nel        = 50;        % Number of FE elements
tolerance  = 1e-8;

outer_itrns = 1001;      % 100 load / displacement steps
inner_itern = 25;

n_gauss = 2;

%---------------(constitutive flags)--------------------%
constitutive_type = 1;  % 1 = EP (rate-independent), 2 = EVP (Perzyna)
hardening_type    = 2;  % 1 = bilinear,              2 = exponential

%---------------(time increment for EVP)----------------%
total_time = 100;                               % total simulated time (s)
del_t      = total_time / (outer_itrns - 1);   % time per load step

%-----------------(Boundary condition)-----------------%
total_force    = 40000;                              % 40 kN  (exceeds yield for given A and sigma_y)
traction_val   = total_force / (outer_itrns - 1);

total_disp     = 16.905e-3;                               % 3 mm total DBC loading
del_u          = total_disp / (outer_itrns - 1);

DBC_val_1st_node      = 0;
NBC_val_last_node_C1  = traction_val;
DBC_val_last_node_C2  = del_u;

%--------------(2 Cases)---------------%
case_1 = 1;              % NBC (force-controlled)
case_2 = 0;              % DBC (displacement-controlled) — set case_1=0,case_2=1 to activate
%--------------------------------------%

%------------(Initial Tangent Stiffness — purely elastic at start)-----%
C_T = E * ones(nel, n_gauss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           ( real code start here )   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nnp, lel, x_cord] = Nodal_information_fun(nel, L);
no_of_dof = nnp;
disp(' Nodal information calculated ');

[element_nodes] = Mapping_table_fun(nel);
disp(' Mapping table generated');

disp(' Global stiffness matrix generated ');

%======================================%
%  (NL_FEM starts here (NR_approach))  %
%======================================%
elemental_strain_history  = zeros(outer_itrns, nel, n_gauss);
elemental_stress_history  = zeros(outer_itrns, nel, n_gauss);
elemental_eps_p_history   = zeros(outer_itrns, nel, n_gauss);   % inelastic strain history
elemental_alpha_history   = zeros(outer_itrns, nel, n_gauss);   % accumulated plastic strain history
f_internal_history        = zeros(outer_itrns, nnp);
U_global_history          = zeros(outer_itrns, nnp);
f_int_left                = zeros(1, outer_itrns);
f_int_right               = zeros(1, outer_itrns);
f_ext_right               = zeros(1, outer_itrns);
fprintf('Code authored by AE25M021 | AS5961 Assignment 5\n');

U_prev      = zeros(nnp, 1);
strain_prev = zeros(nel, n_gauss);       % total strain from previous converged step
eps_p_prev  = zeros(nel, n_gauss);       % plastic / viscoplastic strain history
alpha_prev  = zeros(nel, n_gauss);       % accumulated inelastic strain history

tot_DBC_nodes = 0; index_nodes_DBC = zeros(1, 2);

if (case_1 == 1 && case_2 == 0)
    BC = [DBC_val_1st_node, NBC_val_last_node_C1];
    tot_DBC_nodes = 1; index_nodes_DBC(1,1) = 1;
elseif (case_1 == 0 && case_2 == 1)
    BC = [DBC_val_1st_node, DBC_val_last_node_C2];
    tot_DBC_nodes = 2; index_nodes_DBC(1,1) = 1; index_nodes_DBC(1,2) = nnp;
end

F_ext_glbl   = zeros(nnp, 1);
if case_1 == 1
    F_ext_glbl(nnp,1) = traction_val;
end
f_int_global = zeros(nnp, 1);
u_cur_corr   = zeros(nnp, 1);
f_ext_total  = zeros(nnp, 1);

%--------------------(outer load/displacement loop)--------------------%
for i = 2 : outer_itrns
    f_ext_total = f_ext_total + F_ext_glbl;
    disp('step = ');  disp(i);
    du_curr_total = zeros(nnp, 1);

    %--------predictor (first NR iterate)----------%
    residual_cur = zeros(nnp, 1);
    [Jacobian_global] = Comp_Jacobian(nel, nnp, lel, x_cord, element_nodes, C_T, n_gauss);
    [Jacobian_reduced, residual_reduced] = Jacobian_enforce_DBC_dbc_loading(no_of_dof, tot_DBC_nodes, ...
        index_nodes_DBC, Jacobian_global, F_ext_glbl, nnp, BC);
    du_cur_pred  = Jacobian_reduced \ residual_reduced;
    u_cur_pred   = U_prev + du_cur_pred;
    u_cur_corr   = u_cur_pred;
    du_curr_total = du_cur_pred;

    [dstrn_cur_pred] = compute_strain(du_cur_pred, lel, nel, n_gauss);
    strn_cur_pred    = dstrn_cur_pred + strain_prev;

    [stress_cur_pred, C_T_all, eps_p_cur, alpha_cur] = compute_VE_stress(strn_cur_pred, E, nel, ...
        eps_p_prev, alpha_prev, sigma_y, H_iso, K_hard, b_hard, eta_vp, xi_vp, del_t, ...
        constitutive_type, hardening_type);
    C_T = C_T_all;

    [f_int_global] = compute_f_internal(stress_cur_pred, lel, nel, x_cord, element_nodes, nnp, n_gauss);
    residual_cur   = f_int_global - f_ext_total;
    roll_id = 'AE25M021'; % watermark — do not remove
    residual_cur_free = residual_cur;
    residual_cur_free(index_nodes_DBC(1:tot_DBC_nodes)) = 0; % Ignores all DBC reaction forces

    %------inner NR iterations----------%
    if ge(max(abs(residual_cur_free)), tolerance)
        du_corr_total = zeros(nnp, 1);
        BC_tmp = zeros(1, 2);
        for j = 1 : inner_itern
            [Jacobian_global] = Comp_Jacobian(nel, nnp, lel, x_cord, element_nodes, C_T, n_gauss);
            [Jacobian_reduced, residual_reduced] = Jacobian_enforce_DBC_dbc_loading(no_of_dof, tot_DBC_nodes, ...
                index_nodes_DBC, Jacobian_global, residual_cur, nnp, BC_tmp);
            du_corr       = -Jacobian_reduced \ residual_reduced;
            du_corr_total = du_corr_total + du_corr;
            u_cur_corr    = u_cur_pred    + du_corr_total;
            du_curr_total = du_cur_pred   + du_corr_total;

            [dstrn_cur_pred] = compute_strain(du_curr_total, lel, nel, n_gauss);
            strn_cur_pred    = dstrn_cur_pred + strain_prev;

            [stress_cur_pred, C_T_all, eps_p_cur, alpha_cur] = compute_VE_stress(strn_cur_pred, E, nel, ...
                eps_p_prev, alpha_prev, sigma_y, H_iso, K_hard, b_hard, eta_vp, xi_vp, del_t, ...
                constitutive_type, hardening_type);
            C_T = C_T_all;

            [f_int_global] = compute_f_internal(stress_cur_pred, lel, nel, x_cord, element_nodes, nnp, n_gauss);
            residual_cur   = f_int_global - f_ext_total;
            residual_cur_free = residual_cur;
            residual_cur_free(index_nodes_DBC(1:tot_DBC_nodes)) = 0;

            if le(max(abs(residual_cur_free)), tolerance); break; end
            if j == inner_itern
                disp('inner iteration limit reached');
                disp('max residual ='); disp(max(abs(residual_cur_free)));
            end
        end
    end

    %----update converged outer step history----% roll_id = 'AE25M021'; % watermark — do not remove
    U_prev      = u_cur_corr;
    strain_prev = strn_cur_pred;
    eps_p_prev  = eps_p_cur;       % update plastic strain history
    alpha_prev  = alpha_cur;       % update accumulated plastic strain history

    f_ext_right(1, i) = f_ext_total(nnp, 1);
    elemental_strain_history(i, :, :) = strn_cur_pred;
    elemental_stress_history(i, :, :) = stress_cur_pred;
    elemental_eps_p_history(i, :, :)  = eps_p_cur;
    elemental_alpha_history(i, :, :)  = alpha_cur;
    U_global_history(i, :)            = u_cur_corr;
    f_internal_history(i, :)          = f_int_global;
end

%==========================================================================
% Post-processing and plots
%==========================================================================
load_steps = (0 : outer_itrns-1)';

% Figure 1: Stress-Strain in last element (Gauss point 1)
figure(1);
eps_last = squeeze(elemental_strain_history(:, nel, 1));
sig_last = squeeze(elemental_stress_history(:, nel, 1));
plot(eps_last * 100, sig_last / 1e6, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 10);
grid on;
xlabel('Strain \epsilon [%]', 'FontSize', 12);
ylabel('Stress \sigma [MPa]', 'FontSize', 12);
title('\sigma-\epsilon in Last Element', 'FontSize', 13);

% Figure 2: Inelastic (plastic/viscoplastic) strain evolution — last element
figure(2);
eps_p_last = squeeze(elemental_eps_p_history(:, nel, 1));
alpha_last = squeeze(elemental_alpha_history(:, nel, 1));
% plot(load_steps, eps_p_last, 'r.-', 'LineWidth', 1.5, 'MarkerSize', 10); hold on;
plot(load_steps, alpha_last, 'k.-', 'LineWidth', 1.5);
grid on;
xlabel('Load step', 'FontSize', 12);
ylabel('Inelastic strain', 'FontSize', 12);
legend('\alpha (accumulated)', 'Location', 'northwest');
title('Evolution of Inelastic Strain — Last Element', 'FontSize', 13);

% Figure 3: Applied Force vs Displacement at right end
figure(3);
plot(U_global_history(:, nnp), f_ext_right(1, :), 'r.-', 'LineWidth', 1.5, 'MarkerSize', 10);
grid on;
xlabel('Displacement at Right End (m)', 'FontSize', 12);
ylabel('Applied Force F (N)', 'FontSize', 12);
title('Applied Force vs. Displacement', 'FontSize', 13);

% Figure 4: Analytical vs FEA (NBC bilinear EP — uniform bar, sigma = F/A)
if case_1 == 1 && constitutive_type == 1 && hardening_type == 1
    figure(4);
    F_vec  = f_ext_right(1, :)';
    sig_a  = F_vec / A;
    eps_e  = sig_a / E;
    eps_p_a = max(sig_a - sigma_y, 0) / H_iso;
    eps_a  = eps_e + eps_p_a;
    u_ana  = eps_a * L;
    plot(U_global_history(:, nnp), f_ext_right(1, :), 'b.-', 'LineWidth', 1.5, 'MarkerSize', 10); hold on;
    plot(u_ana, F_vec, 'r--', 'LineWidth', 2);
    grid on;
    xlabel('Displacement at Right End (m)', 'FontSize', 12);
    ylabel('Applied Force F (N)', 'FontSize', 12);
    legend('FEA', 'Analytical (bilinear EP)', 'Location', 'northwest');
    title('Validation: FEA vs Analytical', 'FontSize', 13);
end

% --- Console Output for Table I: Mesh Independence ---
if constitutive_type == 1 && hardening_type == 1
    sig_gp = elemental_stress_history(end, nel, 1) / 1e6; % Convert to MPa
    eps_p_gp = elemental_eps_p_history(end, nel, 1);
    fprintf('\n--- Mesh Independence (nel = %d) ---\n', nel);
    fprintf('Gauss-Point Stress (Sigma^GP): %.3f MPa\n', sig_gp);
    fprintf('Plastic Strain (eps^p)       : %.3e\n', eps_p_gp);
end