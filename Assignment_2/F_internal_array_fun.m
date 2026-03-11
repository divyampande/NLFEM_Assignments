function [ F_internal_global] = F_internal_array_fun (lel, nel, nnp, A, elemental_gauss_point_flux, x_cord, element_nodes  ) 
 [ B] = Shape_function_fun( lel) ;

 F_internal_local = zeros ( 2*nel ,1) ; 
 F_internal_global = zeros(nnp, 1);

%[elemental_strain, elemental_gauss_point_flux]= Gauss_point_flux_fun( nel, lel, U, element_nodes, elemental_pressure, E, beta, del_T, alfa, zeta)
syms  x
for i = 1:nel
    F_internal_local([(2*i-1) , 2*i],:)= int (  B' * A *  elemental_gauss_point_flux(i, 1) , x, x_cord(i), x_cord(i+1) );
end


for i =1:nel
    dof_tmp =  element_nodes(i,:) ;
   % K_global ( element_nodes(i,:), element_nodes(i,:) ) = K_global( element_nodes(i,:), element_nodes(i,:) ) + K_local([(2*i-1) , 2*i], :) ;
    F_internal_global ( dof_tmp, 1 ) =  F_internal_global( dof_tmp ,1 ) + F_internal_local([(2*i-1) , 2*i],:) ;
end
 