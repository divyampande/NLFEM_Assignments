
%========================================================================%
% (1) I m using  (1-point)  gauss quadrature rule for calculating the
% elemental flux.

% (2) Although i have been using (3-point) rule for calculting the
% the integration term involve in the global stiffness matrix
%=========================================================================%

function [elemental_strain, elemental_gauss_point_flux]= Gauss_point_flux_fun( nel, lel, U_curr, element_nodes, elemental_pressure, E, beta, del_T, alfa, zeta)

%disp(nel) ; 
elemental_strain= zeros(nel,1);
elemental_gauss_point_flux= zeros(nel,1);           % Each row stored each element flux value at gauss point

[B] = Shape_function_fun( lel); % Caling derivative of shape function


for i = 1:nel
    elemental_strain(i,1) = B * U_curr(element_nodes(i,:), 1);     % History variable 
    P= elemental_pressure(i,1);                                     % History variable  % each element have one pressure it will go there...
    elemental_gauss_point_flux(i, 1) = E * ( elemental_strain(i,1) -  (beta * del_T) - zeta ) - (alfa * P) ;
end

end

 