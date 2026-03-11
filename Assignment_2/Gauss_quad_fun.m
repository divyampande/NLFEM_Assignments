function [weight_coeff,  gauss_point ]= Gauss_quad_fun(n)
        if  n == 1
            weight_coeff = 2;
            gauss_point = 0;

        elseif  n == 2
            weight_coeff = [1, 1];
            gauss_point = [ 1/sqrt(3), -1/sqrt(3)];

        elseif n == 3
            weight_coeff = [ 5/9,  8/9,  5/9];
            gauss_point= [-sqrt(3/5),  0,  sqrt(3/5)];

        else
            weight_coeff = [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451 ];
            gauss_point = [-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116 ];

        end

    end