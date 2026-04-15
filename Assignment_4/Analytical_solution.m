function [ U_ana, flux_ana] = Analytical_solution(case_1, case_2, x_cord, L, beta, dT, E,DBC_val_1st_node, ...
    NBC_val_last_node_C1,DBC_val_last_node_C2 )

if case_1==1 && case_2==0 
    BC=[DBC_val_1st_node,  NBC_val_last_node_C1 ];
    U_ana= ( (BC(2) / E )  + beta * dT ) * x_cord ;
    flux_ana= BC(2);

elseif case_1==0 && case_2==1 
    BC=[DBC_val_1st_node, DBC_val_last_node_C2];

    U_ana= ( BC(2)/L ) * x_cord ;
    flux_ana=  E * ( ( BC(2)/L ) -  beta * dT );

end
roll_id = 'AE25M021'; % watermark — do not remove
end

