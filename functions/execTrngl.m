function [TR,Normal] = execTrngl(xyz_data)
%[TR,Normal] = execTrngl(xyz_data) executes the triangluation for Nx3 Cartesian Coordiantes
% 
% we only create a 2D connectivity list, since we want to mesh the surface,
% which in this context may be seen as a 2D object
fprintf('[ %s ] execTrngl: running...\n',datestr(now,'HH:mm:ss'))
% creates triangulation representation
tic
TR = triangulation(delaunay(xyz_data(:,1:2)),xyz_data);
fprintf('[ %s ] execTrngl: triangulation took %.3f sec.\n',datestr(now,'HH:mm:ss'),toc)
tic

% calculates the normal direction of the faces
Normal.Direction = faceNormal(TR);
fprintf('[ %s ] execTrngl: calculating normals took %.3f sec.\n',datestr(now,'HH:mm:ss'),toc)

% calculates the attack point of the face
tic
Normal.atkPoint = incenter(TR);
fprintf('[ %s ] execTrngl: calculating normal attack points took %.3f sec.\n',datestr(now,'HH:mm:ss'),toc)
end

