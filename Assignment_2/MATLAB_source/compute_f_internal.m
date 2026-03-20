function [f_int_global] = compute_f_internal(stress_cur_pred, lel, nel, x_cord, element_nodes, nnp, n_gauss)
    [B] = Shape_function_fun(lel);
    F_int_local = zeros(2*nel, 1);
    F_int_glbl = zeros(nnp, 1);
    
    A = 1e-4;
    [weight_coeff, ~] = Gauss_quad_fun(n_gauss);
    detJ = lel / 2;
    
    for i = 1:nel
        f_int_ele = zeros(2, 1);
        
        % Fetch the constant stress for this linear element
        for j = 1:n_gauss
            sigma = stress_cur_pred(i, j);  
            f_int_ele = f_int_ele + weight_coeff(j) * (B' * sigma) * A * detJ;
        end
        
        F_int_local([(2*i-1), 2*i],:) = f_int_ele;
        F_int_glbl(element_nodes(i,:), 1) = F_int_glbl(element_nodes(i,:), 1) + f_int_ele;
    end
    f_int_global = F_int_glbl;
end