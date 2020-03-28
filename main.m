%% Create masked spatial data

stngs = struct();
stngs.plotoutput = true;           % activate plotting of data
stngs.NumFigs = 4;
stngs.ptPlot = false;               % activate calculation and plotting of the point Cloud object
stngs.normalize = true;             % normalize data so a unit of 1, if x/y cell size will be set to 1, z values scaled accordingly
stngs.scalez = 1.5;
stngs.subsamp = true;               % activate subsampling of data
stngs.limitArea = true;
stngs.stlout = false;               % activate exporting of .stl file
stngs.defpath = 'D:\02_Documents\04_Projekte\couchtisch\02_rohdaten\';

%% User specified I/O
importlv = true;                % if no data present, import always
if exist('data','var')          % if data present, let user decide between using it again or loading new data
    switch questdlg('Do you want to import new data?','IMPORT','Yes','No','No')
        case 'No'
            importlv = false;
    end
end
if stngs.stlout
    exportlv = true;                      % initially we always want to catch an export path
    if isfield(stng,'pathout')          % only if an output path exists already, check if user wants to use the same path again
        switch questdlg('Do you want to export to the previous path?','EXPORT','Yes','No','Yes')
            case 'Yes'
                exportlv = false;
        end
    end
    if exportlv
        [filename_out,path_out] = uiputfile('*.stl','Save as:',fullfile(stngs.pin,'*.stl'));
    end
    clearvars export;
end

% Execute Fileimport
if importlv
%     import_data_type = questdlg('What type of data do you want to import?','Filetype','GIS','Matlab','Matlab');
    stngs.dtypes = {'*.mat','MATLAB Files (*.mat)';...
                    '*.asc;*.txt','GIS files (*.asc,*.txt)';...
                    '*.*','All Files (*.*)'};                               % data types
    [stngs.fnin,stngs.pin,stngs.ftind] = uigetfile(stngs.dtypes,'Select data file...',stngs.defpath);
    switch stngs.ftind                                                      % switch behaviour for file type index
        case 0                  % case: abort
            fprintf('[ %s ] user discarded import.\n',datestr(now,'HH:mm:SS'));
            return
        case 1                  % case: Matlab file
            fprintf('[ %s ] loading MATLAB file: %s\n',datestr(now,'HH:mm:SS'),fullfile(stngs.pin,stngs.fnin)); tic;
            load(fullfile(stngs.pin,stngs.fnin));
            fprintf('[ %s ] loading MATLAB file: done. Time: %.3f sec.\n',datestr(now,'HH:mm:SS'),toc);
        case 2                  % case: GIS file
            fprintf('[ %s ] loading GIS file: %s\n',datestr(now,'HH:mm:SS'),fullfile(stngs.pin,stngs.fnin)); tic;
            [data,meta] = arcgridread(fullfile(stngs.pin,stngs.fnin),'planar');
            meta = struct(meta);                                            % convert the arcgis meta struct to a regular struct
            meta.InitProcessDone = false;
            meta.Subsampling.State = false;
            meta.Normalizing.State = false;
            meta.Scale = [1 1 1];                                           % [X Y Z]
            fields2rm = {'DeltaX','DeltaY',...
                        'RasterHeightInWorld','RasterWidthInWorld',...
                        'XWorldLimits','YWorldLimits',...
                        'XIntrinsicLimits','YIntrinsicLimits'};
            meta = rmfield(meta,fields2rm);
            fprintf('[ %s ] loading GIS file: done. Time %.3f sec.\n',datestr(now,'HH:mm:SS'),toc);
    end
end
clearvars importlv;

%% prepare plotting
if stngs.plotoutput
    % create figure
    curFigs = numel(findobj('type','figure'));
    createFigH = false;
    copyFigH = false;
    if ~curFigs                                                             % create the figures "from scratch"
        for fig=1:stngs.NumFigs
            figH(fig) = figure();
        end
    else                                                                    % create our own, known, figure handles, append if neccessary
        figH = findobj('type','figure');
        for fig = (numel(figH)+1):stngs.NumFigs
            figH(fig) = figure();
        end
    end
    for fig=1:stngs.NumFigs
        set(figH(fig),'WindowStyle','docked');
    end
    clearvars fig
end

%% preprocess data after import
if ~meta.InitProcessDone
    data = rmNaN(data);     % remove NaN rows and cols
    meta.RasterSize = size(data);
    meta = updateGISMetadata(meta);
    meta.InitProcessDone = true;
end
if and(stngs.normalize,~meta.Normalizing.State)
    if meta.CellExtentInWorldX == meta.CellExtentInWorldY           % check for uniformity of original data scaling
        meta.Normalizing.Factor = meta.CellExtentInWorldX;
    else
        % we need to determine which axis we want everything to reference to
    end
    data = data - min(data,[],'all');                               % remove z-axis offset
    data = data / meta.Normalizing.Factor;                         % to keep physical reference, we scale the values by the scaling factor of the grid
    meta.CellExtentInWorldX = meta.CellExtentInWorldX / meta.Normalizing.Factor;
    meta.CellExtentInWorldY = meta.CellExtentInWorldY / meta.Normalizing.Factor;
    meta.Transformation.Jacobian.Numerator = meta.Transformation.Jacobian.Numerator / meta.Normalizing.Factor;
    meta = updateGISMetadata(meta);
    meta.Normalizing.State = true;
end
if ~(stngs.scalez / meta.Scale(3) == 1)                             % scale data in z to match the set scale factor
    fprintf('max value before scaling: %.2f\n',max(max(data)))
    data = data * (stngs.scalez / meta.Scale(3));
    meta.Scale(3) = meta.Scale(3)  * (stngs.scalez / meta.Scale(3));
    fprintf('max value after scaling: %.2f\n',max(max(data)))
end
if and(stngs.subsamp,~meta.Subsampling.State)                       % reduce data, only if requested and not previously done
%     tic;
%     data = downsample(downsample(data',3)',3);                              % subsampling
%     fprintf('[ %s ] subsampling took: %.3f sec.\n',datestr(now,'HH:mm:ss'),toc);
    [data, meta] = execSubsample(data,meta,3);
end


% Create grid
x_vals =  (0:size(data,2)-1)' * meta.CellExtentInWorldX;
y_vals =  flipud((0:size(data,1)-1)') * meta.CellExtentInWorldY;

%% limit processed area
xmask = ones(1,size(data,2));   % row vector
ymask = ones(size(data,1),1);   % col vector
if stngs.limitArea
    xmask = xmask * 0;
    ymask = ymask * 0;
    areaname = 'untersee_subsection';
    switch areaname
        case 'surface'
            x_lim = 3950:4050;
            y_lim = 1950:2050;
        case 'untersee'
            x_lim = 1:1000;
            y_lim = 1400:2100;
        case 'untersee_subsection'
            x_lim = 1:200;
            y_lim = 1650:2100;
    end
    if meta.Subsampling.State
        x_lim = ceil(x_lim / meta.Subsampling.Factor);
        y_lim = ceil(y_lim / meta.Subsampling.Factor);
    end
    %     data = data(y_lim,x_lim);
    %     x_vals = x_vals(x_lim);
    %     y_vals = y_vals(y_lim);
    xmask(x_lim) = 1;
    ymask(y_lim) = 1;
    ymask = ymask;
end
xmask = logical(xmask);
ymask = logical(ymask);
datamask = logical(xmask .* ymask);
xyz_data=grid2cart(x_vals(xmask),y_vals(ymask),data(datamask));

set(0,'CurrentFigure',figH(1))
figH(1).Name = 'ImshowPair';
imshowpair(data,datamask)

clearvars areaname datamask x_lim y_lim xmask ymask

return;
%% Triangulation
[TR,Normal] = execTrngl(xyz_data);

%% using alpha shape
TR2 = execAlphaShp(xyz_data,'PlotOut',true,'HoleThreshold',1e4);
umfang = TR2.freeBoundary;
% find distance distribution
[~,distvec] = nearestNeighbor(TR2,TR2.Points);
pcshow(TR2.Points,squeeze(ind2rgb((distvec./max(distvec,[],'all')),parula)));

if stngs.ptPlot
    tic
    ptCloud = pointCloud(xyz_data);                                         % it might not be neccessary to create a point cloud fist
    ptCloud.Normal = pcnormals(ptCloud);                                    % ToDo: performance tests
    fprintf('[ %s ] create ptCloud object: %.3f sec.\n',datestr(now,'HH:mm:ss'),toc)
end
%% plotting
if stngs.plotoutput
    set(0,'CurrentFigure',figH);        % activate current figure
    clf; hold on;
    % point cloud of raw data
    if stngs.ptPlot
        pcshow(ptCloud);
    end
    % normal vectors of surfaces
    quiver3(Normal.atkPoint(:,1),Normal.atkPoint(:,2),Normal.atkPoint(:,3),Normal.Direction(:,1),Normal.Direction(:,2),Normal.Direction(:,3));
    % attack points of normal vectors
    scatter3(Normal.atkPoint(:,1),Normal.atkPoint(:,2),Normal.atkPoint(:,3),'.r');
    % plot the mesh
    trimesh(TR)
    % axes setup
    if ~stngs.limitArea
        x_lim = [3950, 4050];
        y_lim = [1950, 2050];
        if meta.Subsampling.State
            x_lim = ceil(x_lim / meta.Subsampling.Factor);
            y_lim = ceil(y_lim / meta.Subsampling.Factor);
        end
        set(gca,'XLim',x_vals(x_lim),'YLim',flipud(y_vals(y_lim)));
    end
end

%% make .stl
if stngs.stlout
    fname.inbuild = ['inbuild_' filename_out];
    fname.fex = ['fex_' filename_out];
    
    tic;
    stlwriteFEX(fullfile(path,fname.fex),TR.ConnectivityList,TR.Points);
    fprintf('[ %s ] stlwriteFEX took: %.3f sec.\n',datestr(now,'HH:mm:ss'),toc)
    
    tic;
    stlwrite(TR,fullfile(path,fname.inbuild))
    fprintf('[ %s ] inbuild stlwrite took: %.3f sec.\n',datestr(now,'HH:mm:ss'),toc)
    % ToDo:
    % test speed against the new inbuild function
end
% clear x_vals...
%     y_vals;