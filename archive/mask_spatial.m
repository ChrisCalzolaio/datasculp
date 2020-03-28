function [mask] = mask_spatial(spatial_data,no_data,disp_image,disp_log)
% MASK_SPATIAL Mask spatial data according to no-data value
%
%   MASK_SPATIAL(spatial_data,-9999,1,1) masks all cells with value -9999
%   with 0, plots the mask and prints the logs
%
%   mask values:
%   0:  empty data cell
%   1:  payload cell
%   2:  no-data cells that you don't want masked (neighbours, border of
%   data

[nrows,ncols] = size(spatial_data);
mask = single(zeros(nrows,ncols));
npayload = 0;
nneighbour = 0;
nborder = 0;
nwhite = 0;

% progressbar
num_calc = nrows * ncols;
h = waitbar(0,sprintf('Operation %d of %d',0,num_calc));
n_waitbar = 1;

for l=1:nrows
   for k=1:ncols
       if spatial_data(l,k) > no_data
           mask(l,k) = 1;
           npayload = npayload + 1;
       elseif l == 1 || l == nrows || k == 1 || k == ncols
           mask(l,k) = 2;
           nborder = nborder + 1;
       else
           [box] = spatial_data(l-1:l+1,k-1:k+1);
           check = box == no_data;
           check = all(all(check));
           if check == 0
               mask(l,k) = 2;
               nneighbour = nneighbour + 1;
           else
               mask(l,k) = 0;
               nwhite = nwhite + 1;
           end
       end
   end
   
   % improve runtime 
   if n_waitbar == 1000
       waitbar((l*k)/num_calc,h,sprintf('Operation %d of %d',l*k,num_calc));
       n_waitbar = 1;
   end
   n_waitbar = n_waitbar + 1;
   
end



if disp_image
    image(mask,'CDataMapping','scaled')
%    image(mask)
end
if disp_log
    fprintf('Number of payload data: %d\nNumber of neigbours: %d\nNumber of border points: %d\nNumber of empty cells: %d\n',npayload,nneighbour,nborder,nwhite)
end

close(h); clear h;

end