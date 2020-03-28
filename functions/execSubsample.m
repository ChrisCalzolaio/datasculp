function [data,meta] = execSubsample(data,meta,n)
% [data,meta] = execSubsample(data,meta,n) subsamples data and updates
% metadata
tic; fprintf('[ %s ] execSubsample: running...\n',datestr(now,'HH:mm:ss'))
data = downsample(downsample(data',n)',n);                                  % subsampling

meta.Subsampling.Factor = n;
meta.Subsampling.OrigRasterSize = meta.RasterSize;
meta.Subsampling.OrigCellExtentInWorldX = meta.CellExtentInWorldX;
meta.Subsampling.OrigCellExtentInWorldY = meta.CellExtentInWorldY;

meta.CellExtentInWorldX = meta.CellExtentInWorldX * n;                      % when downsampling, cell size effectively is increased
meta.CellExtentInWorldY = meta.CellExtentInWorldY * n;
meta.RasterSize = ceil(meta.RasterSize ./ n);
meta.Transformation.Jacobian.Numerator = meta.Transformation.Jacobian.Numerator * n;

meta = updateGISMetadata(meta);

fprintf('[ %s ] execSubsample: subsampling took: %.3f sec.\n',datestr(now,'HH:mm:ss'),toc);
meta.Subsampling.State = true;
end