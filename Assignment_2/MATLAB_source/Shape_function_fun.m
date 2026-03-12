function [ B] = Shape_function_fun( lel)

  
  
    %N=[ ((x_b-x)/(x_b - x_a)) , ((x-x_a)/(x_b - x_a))   ]; % defined linear shape function-
    B= [-1/lel, 1/lel];
    
end 