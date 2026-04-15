% damage is dependenet on disp so no need to update in main file-
% the moment we update the disp damage will automatically will get updated

%------------------------------------------------------
% This function is just to update stress without updating damage
%-------------------------------------------------------
% function [elmn_sig, elmn_strn]= stress_update_fun(U_prev, lel, nel, E, D_cur)
% [B] = Shape_function_fun( lel);
% elmn_strn= zeros(nel,1);  elmn_sig= zeros(nel,1);
%
% %---------------------(obtaining elemental stress)----------------------%
% for i = 1:nel
%     elmn_strn(i, 1)= B*U_prev([i,i+1],1);                       % elemental strain
%     elmn_sig(i, 1)= (1- D_cur(i, 1))*E*elmn_strn(i, 1);     % elemental stress
% end

function [elmn_sig, elmn_strn, D_cur]= stress_update_fun(U_prev, lel, nel, W, E, m, D_tr)
[B] = Shape_function_fun( lel);
elmn_strn= zeros(nel,1); elmn_strn_enrgy= zeros(nel,1); elmn_sig= zeros(nel,1);

%---------------------(obtaining elemental stress)----------------------%
for i = 1:nel
    elmn_strn(i, 1)= B*U_prev([i,i+1],1);                       % elemental strain
    elmn_strn_enrgy(i, 1)= (1/2)*E*elmn_strn(i)^2;              % elemental strain energy
    f_tr(i, 1)=(1-D_tr(i, 1))^(1/ m)* elmn_strn_enrgy(i, 1)- W; % surface equation
    %------------------------------------------------
    if (le(f_tr(i, 1),0))                                       % inside yield surface
        disp('------');
        D_cur(i, 1) = D_tr(i, 1);                               % (update damage)
    else                                                        % on the surface f_tr==0
        disp('+++++++')
        D_cur(i, 1)= 1- (W/elmn_strn_enrgy(i, 1))^(1/m);        % (update damage)
    end
    elmn_sig(i, 1)= (1- D_cur(i, 1))*E*elmn_strn(i, 1);     % elemental stress
    %--------------------------------------------------
end
