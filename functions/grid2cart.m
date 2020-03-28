function [cartCoord] = grid2cart(xvec,yvec,hmap)
%[cartCoord] = grid2cart(xvec,yvec,hmap) Creates a Nx3 matrix of cartesian coordinates from gridded data
%   expects empty cell values to be replaced by NaNs beforehand

tic
xvec = xvec(:)';                        % returns a row vector
yvec = yvec(:);                         % returns a col vector
msk = ~isnan(hmap);                    % masks NaN data out, true if data point is payload data, false if NaN
[xmsh,ymsh] = meshgrid(xvec,yvec);  % we create a mesh of X and Y vector, dimensions according to the heightmap
cartCoord = [xmsh(msk(:)),ymsh(msk(:)),hmap(msk(:))];    % we create a columnwise vector of the X, Y and Z data
fprintf('[ %s ] grid2cart: run time was: %.3f sec.\n',datestr(now,'HH:mm:ss'),toc)
end

