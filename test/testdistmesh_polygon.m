datasource = 'bodensee';
hulltype = 'alpha';

switch datasource
    case 'bodensee'
        x = xyz_data(:,1);
        y = xyz_data(:,2);
    case 'rand'
        x = rand(1,1000)';
        y = rand(1,1000)';
        scatter(x,y)
    case 'example'
        pv=[-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
        x = pv(:,1);
        y = pv(:,2);
end

switch hulltype
    case 'alpha'
        [TR,shp] = execAlphaShp([x,y],'Holethreshold',1e4);
        k = TR.freeBoundary;
    case 'conv'
        k = convhull(x,y);
end
umfangx = x(k);
umfangy = y(k);
pv = [umfangx(:), umfangy(:)];
dx = diff(umfangx(:));
dy = diff(umfangy(:));
length = sqrt(dx.^2 + dy.^2);
bbox = [floor(min(pv)); ceil(max(pv))];
bbox = [-100, 1850;300, 2100];

numFigs = 4;
if ~exist('figH','var')
    for figs = 1:numFigs
        figH(figs) = figure();
    end
elseif string(class(figH)) == 'matlab.ui.Figure'
    for figs = 1:min(numel(figH),numFigs)
        if ~isvalid(figH(figs))
            figH(figs) = figure();
        end
    end
    for figs=1:(numFigs-numel(figH))
        figH(numel(figH)+1) = figure();
    end
end

[p,t]=distmesh2dDiscrete(shp,min(length),bbox,figH);