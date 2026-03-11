function [nnp, lel, x_cord] = Nodal_information_fun(nel, L)
nnp= nel+1;            % number of nodal point
lel= L/nel;            % length of each element( having 3 node)
x_cord= 0:lel:L;       % x co-ordinate
end