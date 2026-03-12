function [K_global]= Stiffness_matrix_fun( nel, nnp, lel,  A, E,  x_cord,  element_nodes)

K_local =zeros (2*nel,2) ;
K_global = zeros(nnp, nnp);

[B] = Shape_function_fun( lel);     % Caling derivative of shape function






for i = 1:nel
    syms x
    K_local([(2*i-1) , 2*i], :) = int ( (A*E)*( B' * B ), x, x_cord(i), x_cord(i+1));     % local stiffness in each loop
end

%--------------(Assembly of global stiffness matrix)-----------------%
for i =1:nel
    K_global ( element_nodes(i,:), element_nodes(i,:) ) = K_global( element_nodes(i,:), element_nodes(i,:) ) + K_local([(2*i-1) , 2*i], :) ;
end


end