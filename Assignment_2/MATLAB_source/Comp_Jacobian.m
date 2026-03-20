%function to compute Jacobian matrix for VE case
function [Jacobian_global] = Comp_Jacobian(nel, nnp, lel, x_cord, element_nodes, C_T)
    Jacobian_local = zeros(2*nel, 2);
    Jacobian_global = zeros(nnp, nnp); 
    [B] = Shape_function_fun(lel);     % Calling derivative of shape function
    
    A = 1e-4; % Cross-sectional area
    n_gauss = 2; % 2-point Gauss quadrature
    [weight_coeff, ~] = Gauss_quad_fun(n_gauss);
    detJ = lel / 2; % Jacobian of mapping for 1D element
    
    for i = 1:nel
        K_local = zeros(2, 2);
        dof_tmp = element_nodes(i, :);
        for j = 1:n_gauss
            E_t = C_T(i, j);              
            K_local = K_local + weight_coeff(j) * (B' * E_t * B) * A * detJ;
        end
        
        Jacobian_local([(2*i-1), 2*i],:) = K_local;
        Jacobian_global(dof_tmp, dof_tmp) = Jacobian_global(dof_tmp, dof_tmp) + K_local;
    end
end