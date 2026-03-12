function [conv_elmn_strn, conv_elmn_sig]=stress_strain_conv(lel,nel, U_prev, E, W, m)

[B] = Shape_function_fun( lel);                          
conv_elmn_strn= zeros(nel,1);conv_elmn_sig= zeros(nel,1); 

for i = 1:nel
    conv_elmn_strn(i, 1)= B*U_prev([i,i+1],1);              % elemental strain  
    elmn_strn_enrgy(i, 1)= (1/2)*E*conv_elmn_strn(i)^2;     % elemental strain energy
    if (elmn_strn_enrgy(i,1)==0)                            % setting for denumr zero
        elmn_strn_enrgy(i, 1)= 1;                           % setting for denumr zero                          
    end
    elmn_damage(i, 1)= 1-(W/elmn_strn_enrgy(i))^(1/m);      % elemental damage
    conv_elmn_sig(i, 1)= (1- elmn_damage(i, 1))*E*conv_elmn_strn(i, 1); % elemental stress
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%