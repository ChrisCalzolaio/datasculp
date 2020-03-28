function [spatial_data_masked] = apply_mask(spatial_data,mask,repl_value)
%APPLY_MASK Create clean spatial data from raw data and mask

% APPLY_MASK(data,mask,400) replaces all cell masked 0 with NaN, all cells
% masked 2 with the replacement value, all masked 1 are left as is

[nrows,ncols] = size(spatial_data);

spatial_data_masked = nan(nrows,ncols);

for l = 1:nrows
    for m = 1:ncols
        if mask(l,m) == 1
            spatial_data_masked(l,m) = spatial_data (l,m);
        elseif mask(l,m) == 2
            spatial_data_masked(l,m) = repl_value;
        end
    end
end