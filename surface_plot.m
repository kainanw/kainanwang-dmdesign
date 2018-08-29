function [PV,RMS]=surface_plot(x,y,z,diameter,plot_option,remove_ppt,view_angle)

x_grid=linspace(-diameter/2,diameter/2,100);
y_grid=linspace(-diameter/2,diameter/2,100);

[X_grid,Y_grid]=meshgrid(x_grid,y_grid);

X_grid_column=reshape(X_grid,size(X_grid,1)*size(X_grid,2),1);
Y_grid_column=reshape(Y_grid,size(Y_grid,1)*size(Y_grid,2),1);

%% Not_NaN_index define

mask_not_nan=ones(size(X_grid_column));

for ii=1:length(X_grid_column)
    distance=sqrt(X_grid_column(ii)^2+Y_grid_column(ii)^2);
    if(distance>diameter/2*0.95)
       mask_not_nan(ii)=NaN;
    end
    clear distance
end
Grid_mask=reshape(mask_not_nan,size(X_grid,1),size(X_grid,2));

% X_grid_column_not_nan=X_grid_column(isfinite(mask_not_nan));
% Y_grid_column_not_nan=Y_grid_column(isfinite(mask_not_nan));

X_grid=X_grid.*Grid_mask;
Y_grid=Y_grid.*Grid_mask;

%% Intepolation

F=scatteredInterpolant(x,y,z,'natural','nearest');
Z_grid=F(X_grid,Y_grid);

if(remove_ppt==1)
    [coefs]=LS_fit_plan_v2(X_grid,Y_grid,Z_grid);	% Fitting of the best plan
    W_PTT=coefs(1)+coefs(2)*X_grid+coefs(3)*Y_grid;	% Best plan: Surface
    Z_grid=Z_grid-W_PTT;
    clear W_PTT;
end

PV=max(max(Z_grid))-min(min(Z_grid));
RMS=get_RMS(Z_grid,[]);

if(plot_option==1)
    surf(X_grid,Y_grid,Z_grid);
    view(view_angle(1),view_angle(2));
    shading interp
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('W [m]');

    xlim([-diameter/2,diameter/2]);
    ylim([-diameter/2,diameter/2]);
    axis square;
    colormap(jet);
%     colorbar;
end
