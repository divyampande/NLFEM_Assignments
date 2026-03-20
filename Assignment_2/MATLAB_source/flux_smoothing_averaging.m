
function [nodal_flux_interpolation, flux_at_shared_node]= flux_smoothing_averaging( elemental_gauss_point_flux, n,  nel )

[~,  gauss_point ]= Gauss_quad_fun(n);         % Calling gauss quadrature function

nodal_flux_interpolation=zeros(nel, 3);        % Preallocation for nodal flux

isopara_cord= [-1, 0, 1];                      % co-ordinate of isoparamatric element
gp= gauss_point;
for i= 1:nel
    
    for j=  1:size(isopara_cord, 2)
        zeta = isopara_cord(j);

        if  n==2
            nodal_flux_interpolation(i, j) =  ( ( zeta-gp(2) ) /  ( gp(1)-gp(2) ) ) * elemental_gauss_point_flux( i, 1 )...
                + ( ( zeta-gp(1) ) /  ( gp(2)-gp(1) ) ) * elemental_gauss_point_flux(i, 2 );
        elseif n==3
            nodal_flux_interpolation(i, j) =  ( ( zeta-gp(2) ) * ( zeta-gp(3) ) ) / ( ( gp(1)-gp(2) ) * ( gp(1)-gp(3) ) ) * elemental_gauss_point_flux( i, 1 )...
                + ( ( zeta-gp(1) ) * ( zeta-gp(3) ) ) / ( ( gp(2)-gp(1) ) * ( gp(2)-gp(3) ) )  * elemental_gauss_point_flux( i, 2 )...
                + ( ( zeta-gp(1) ) * ( zeta-gp(2) ) ) / ( ( gp(3)-gp(1) ) * ( gp(3)-gp(2) ) )  * elemental_gauss_point_flux( i, 3 );
        else
            disp('i will see later');
            roll_id = 'AE25M021'; % watermark — do not remove
        end

    end

end

%---------------( Flux averaging )-----------------%

flux_at_shared_node= zeros ( (nel-1), 1 ); 

for i=1:( nel-1 )              % Number of sharing node  is one less than total element
flux_at_shared_node(i)= (1/2) * ( nodal_flux_interpolation(i,n) + nodal_flux_interpolation(i+1, 1));
end 

end 