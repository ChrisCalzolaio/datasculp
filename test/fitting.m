tic
import = 0;

data_path = 'D:\02_Documents\04_Projekte\couchtisch\matlab\daten';
data_file = '10m.asc';
data_name = fullfile(data_path,data_file);

if import
    [data,meta] = arcgridread(data_name,'planar');
    x_vals =  (0:meta.RasterSize(2)-1) * meta.CellExtentInWorldX;
    y_vals =  (0:meta.RasterSize(1)-1)' * meta.CellExtentInWorldY;
end

% bounds for 2D testing
run_2d = 0;
if run_2d
    row = 2000;
    x_start = 2699;
    x_end = 5045;
    x_sel = x_vals(x_start:x_end);
    y_sel = data(row,x_start:x_end);
    % show selected data
    plot(x_sel,y_sel,'.-')
end
% bounds for 3D testing
run_3d = 1;
if run_3d
    x_start = 3650;
    x_end = 4100;
    y_start = 1900;
    y_end = 2100;
    
    
    x_sel = x_vals(x_start:x_end);
    y_sel = y_vals(y_start:y_end);
    z_sel = data(y_start:y_end,x_start:x_end);
    surf(x_vals,y_vals,data*15)
    xlabel('x');ylabel('y');zlabel('z');
    shading interp
%     [x_grid,y_grid] = meshgrid(x_vals,y_vals);
%     plot3(x_grid,y_grid,z_sel,'.')
    xlabel('x');ylabel('y');zlabel('z');
    grid; grid minor;
    view([0 0 180])
end

%CREATEFIT(X_VALS,Y_VALS,DATA)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : x_vals
%      Y Input : y_vals
%      Z Output: data
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 07-Apr-2019 21:06:11


%% Fit: 'untitled fit 1'.
[xData, yData, zData] = prepareSurfaceData( x_vals, y_vals, data );

% Set up fittype and options.
ft = 'nearestinterp';

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult);
legend( h, 'untitled fit 1', 'data vs. x_vals, y_vals', 'Location', 'NorthEast' );
% Label axes
xlabel x_vals
ylabel y_vals
zlabel data
grid on
view( -222.1, 26.6 );

toc