function [ P_cur ] = NR_fun_1D ( alfa, phi_1, E, sigma, beta, del_T, lambda, P_tr )

syms P         % f_1 will be function of 'P' 

%------------function in term of syms and subfunction-----------%
pvs=(alfa/(phi_1*E)) * sigma + ( (1+phi_1)/(phi_1)) * ((alfa*P)/E ) - (P/E) + beta*del_T;
f_1 = P*(1+pvs) - lambda ;


jac = jacobian(f_1, P);     % jacobian in term of syms

%---------Initialisation of syms--------------%
P=P_tr;     % calc jacob at P_tr for first time          % for next loop P is not anymore syms , so we need to define syms inside the loop.
X = P ;     % for the first iteration X = P_tr

%-------------while loop----------------------%
tolerance = 1e-12;
iter = 100;
z = 1;
error = 1;

while(ge(error, tolerance) == 1)            % residual calculation-
    Res = eval(f_1);
    jacobian_val = eval(jac);               % jacobian calculation-
    
    
    X_cur = X - ( 1/(jacobian_val) )*Res;   % current predicted value array-
    diff = X_cur - X;                       % compute diff-

    %nr = 0; dr = 0;                        % for single variable we dont need this-
    nr =diff^2;
    dr = X^2;                               % Last predicted value 

    if (dr==0)   % we are setting  the de-numerator to be non zero, as we don't want any undefined situation-
        dr =1;
    end

    error = sqrt(nr/dr);                  % error calculation-
    X = X_cur ;                           % update current vals with predicted vals-
    P=X ;                                 % extract syms from updated current vals
    z = z + 1;                            % cheak for iteration-
    if z > iter
        disp ('Convergence not achieved');
        break;
    end
end
P_cur = P ; 

end


