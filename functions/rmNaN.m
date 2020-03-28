function [data] = rmNaN(data)
%[data] = rmNaN(data) Removes blank rows and columns
%
%   blank meaning only populated with NaN

tic;
compl = isnan(data);                        % complete mask of NaN values
colrow = or(all(compl,2),all(compl,1));     % masks out columns and rows entirely made out of NaN
croppeddim = [numel(find(~all(compl,2))),numel(find(~all(compl,1)))];   % find dimension of cleansed griddata
data = reshape(data(~colrow),croppeddim);   % recreate the matrix structure of the data

fprintf('[ %s ] rmNaN: run time was: %.3f sec.\n',datestr(now,'HH:mm:ss'),toc)
end

