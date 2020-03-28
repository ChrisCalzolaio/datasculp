function make_stl(tri_LUT,xyz_data,path)
%MAKE_STL Writes stl file using a triangulation matrix of m x 3

fid = fopen(path,'wb+');

% count facets
nfacets = size(tri_LUT,1);

% progress bar
h = waitbar(0,sprintf('Facet %d of %d',0,nfacets));
n_waitbar = 0;

% write stl file
title_str = sprintf('Created by make_stl.m %s',datestr(now));
str = sprintf('%-80s',title_str);    
fwrite(fid,str,'uchar');         % Title
fwrite(fid,nfacets,'int32');     % Number of facets

stl_char = nan(nfacets*5,3);

for l=1:nfacets
    p1 = xyz_data(tri_LUT(l,1),:);
    p2 = xyz_data(tri_LUT(l,2),:);
    p3 = xyz_data(tri_LUT(l,3),:);
    local_write_facet(fid,p1,p2,p3);
    
    if n_waitbar == 1000
       waitbar(l/nfacets,h,sprintf('Operation %d of %d',l,nfacets));
       n_waitbar = 0;
   end
   n_waitbar = n_waitbar + 1;
end


fclose(fid);
fprintf('Wrote %d facets',nfacets);
close(h); clear h;

% Local subfunctions
% partly stolen from: surf2stl by Bill McDonald
function local_write_facet(fid,p1,p2,p3)

    v1 = p2-p1;
    v2 = p3-p1;
    v3 = cross(v1,v2);
    n = v3 ./ sqrt(sum(v3.*v3));
    
    fwrite(fid,n ,'float32');
    fwrite(fid,p1,'float32');
    fwrite(fid,p2,'float32');
    fwrite(fid,p3,'float32');
    fwrite(fid,0 ,'int16'  );  % unused
end
end