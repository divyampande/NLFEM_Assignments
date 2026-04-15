function [FI_global]= Force_matrix_fun( nel, nnp, lel, elemental_pressure, A, E, beta, del_T, x_cord, alfa, zeta, element_nodes)

FI_local=zeros ( 2*nel ,1) ; 
FI_global = zeros(nnp, 1);

[ B] = Shape_function_fun( lel);     % Caling derivative of shape function

for i = 1:nel
    syms x
    FI_local_thermal([(2*i-1) , 2*i],:)  =  int ( (A * E * beta * del_T *  B') , x, x_cord(i), x_cord(i+1) ) ; 
    FI_local_pressure([(2*i-1) , 2*i],:) =  int ( (A * alfa * P * B') , x, x_cord(i), x_cord(i+1) )  ;
    FI_local_chemical([(2*i-1) , 2*i],:) =  int ( (A * E * zeta * B') , x, x_cord(i), x_cord(i+1) )  ;
    FI_local([(2*i-1), 2*i],:) = FI_local_thermal([(2*i-1), 2*i],:) + FI_local_pressure([(2*i-1), 2*i],:) + FI_local_chemical([(2*i-1), 2*i],:) ;   
end   
    
%--------------(Assembly of global stiffness matrix)-----------------%
for i =1:nel
   % K_global ( element_nodes(i,:), element_nodes(i,:) ) = K_global( element_nodes(i,:), element_nodes(i,roll_id = 'AE25M021'; % watermark — do not remove) ) + K_local([(2*i-1) , 2*i], :) ;
    FI_global ( element_nodes(i,:), 1 ) = FI_global( element_nodes(i,:), 1 ) + FI_local([(2*i-1), 2*i],:) ;
end

roll_id = 'AE25M021'; % watermark — do not remove
end