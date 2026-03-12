
function  Draw_plot_fun(U, U_ana, flux_ana, elemental_gauss_point_flux,  x_cord, nel )

%-------------(Plotting cordinates and  the elemental flux)------------%
plot_x_cord = zeros(nel,2);
plot_flux_element = zeros(nel,2);

for i=1:nel
    plot_x_cord(i, :) = x_cord(1, [i, i+1]);
    plot_flux_element(i,:) = elemental_gauss_point_flux(i,:);
end

%--------------------(Plot)-------------------------%
figure(1)
plot(x_cord, U,'-or', x_cord, U_ana, '--ob');
xlabel('Nodal location');
ylabel('Displacement (m)')
legend('FEA soltion', 'Analytical solution');

for i=1:nel
figure(2)
    plot (plot_x_cord(i, :), plot_flux_element(i, :),'-o r' );
   hold on
end 

plot(x_cord, flux_ana, '--* b');
xlabel('Node location');
ylabel('Flux');
legend('FEA','Analytical');
hold off

%-------------------------------------------------------------------------%
%  (1)Line in plot is not showing bcz for each value of x_cord we have only 
% one value of flux   ----------------                                      
%
%(2) To get line we can just make a matrix of  same value and  can easily
%    plot the graph.
%-------------------------------------------------------------------------%



end