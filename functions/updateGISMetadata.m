function meta = updateGISMetadata(meta)
% meta = updateGISMetadata(meta)
% update dependent metadata entries

meta.Intrinsic.NumRows = meta.RasterSize(1);
meta.Intrinsic.NumColumns = meta.RasterSize(2);

meta.RasterExtentInWorldY = meta.CellExtentInWorldY * meta.RasterSize(1);
meta.RasterExtentInWorldX = meta.CellExtentInWorldX * meta.RasterSize(2);

meta.YLimWorld = meta.Transformation.TiePointWorld(2) + [0 meta.RasterExtentInWorldX];
meta.XLimWorld = meta.Transformation.TiePointWorld(1) + [0 meta.RasterExtentInWorldX];

meta.YLimIntrinsic = meta.Transformation.TiePointIntrinsic(1) + [0 meta.RasterSize(1)];
meta.XLimIntrinsic = meta.Transformation.TiePointIntrinsic(2) + [0 meta.RasterSize(2)];

end