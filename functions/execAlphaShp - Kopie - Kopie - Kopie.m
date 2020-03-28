function [TR,shp] = execAlphaShp(p,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%
%   P: points, Nx2 or Nx3 matrix of vertices
%
fprintf('[ %s ] execAlphaShp: running...\n',datestr(now,'HH:mm:ss'))
oat = tic;
mode = size(p,2);                % 2d or 3d mode
lvPlot = 0;                 % logic value Plotting
    alphaShpArgs = {};
if nargin>1
    for ind = 1:2:numel(varargin)
        switch lower(varargin{:,ind})
            case 'plotout'
                switch lower(varargin{ind+1})
                    case 'all'
                        lvPlot = 3;
                    case '2d'
                        lvPlot = 2;
                    case '3d'
                        lvPlot = 1;
                end
            case 'holethreshold'
                alphaShpArgs = {'HoleThreshold',varargin{ind+1}};
            otherwise
                error('unrecognized Command');
        end
    end
end

x = p(:,1);                             % single col vector of x Coordinates
y = p(:,2);                             % single col vector of y Coordinates
tic
shp = alphaShape(x,y,alphaShpArgs{:});                  % create the concave Shape around the data points
fprintf('[ %s ] execAlphaShp: alphaShape took %.3f sec.\n',datestr(now,'HH:mm:ss'),toc);tic
alphaTR = alphaTriangulation(shp);      % create a connectivity list of the faces within the shape
fprintf('[ %s ] execAlphaShp: alphaTriangulation took %.3f sec.\n',datestr(now,'HH:mm:ss'),toc);tic
TR = triangulation(alphaTR,p);          % create the triangulation object for ouput
fprintf('[ %s ] execAlphaShp: triangulation took %.3f sec.\n',datestr(now,'HH:mm:ss'),toc);
if lvPlot                               % the user may specify, wether or not to get plots
    cla;
    if lvPlot>1                         % plot 2D only
        plot(shp)
        shpplt.xlim = get(gca,'XLim');
        shpplt.ylim = get(gca,'YLim');
        hold on;
        plot(x,y,'.r')
    end
    if any([1 3] == lvPlot)
        switch mode
            case 2
                trimesh(TR.ConnectivityList,TR.Points(:,1),TR.Points(:,2))
            case 3
                trimesh(TR)
                view(3)
        end
    end
    if lvPlot == 3
        set(gca,'Xlim',shpplt.xlim,'YLim',shpplt.ylim,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1]);
    end
    hold off
end
fprintf('[ %s ] execAlphaShp took %.3f sec.\n',datestr(now,'HH:mm:ss'),toc(oat));
end