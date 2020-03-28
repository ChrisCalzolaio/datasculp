pv=[-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
shp = alphaShape(pv(:,1),pv(:,2),2.5);
alphaTR = alphaTriangulation(shp);      % create a connectivity list of the faces within the shape
TR = triangulation(alphaTR,pv);          % create the triangulation object for ouput
plot(pv(:,1),pv(:,2),'.r');
hold on
plot(shp);

h0 = 0.1;

[x,y]=meshgrid(floor(min(pv(:,1))):h0:ceil(max(pv(:,1))),...
    floor(min(pv(:,2))):h0/sqrt(3)/2:(ceil(max(pv(:,2)))+.5));
x = x(:);y=y(:);

pIN = inShape(shp,x,y);

p = [x(~pIN),y(~pIN)];

plot(p(:,1),p(:,2),'.b')

ind = nearestNeighbor(shp,p(:,1),p(:,2));

plot(p(1,1),p(1,2),'*r')
plot(shp.Points(ind(1),1),shp.Points(ind(1),2),'*r')

con = shp.Points(ind,:) - p;            % connections
conl = vecnorm(con,2,2);                % connection lengths
cond = con./conl;                       % connection directions
q = quiver(p(:,1),p(:,2),cond(:,1),cond(:,2),0.5);

pkeep = ~(rand(size(conl,1),1)<conl);   % points to keep
p = p(pkeep,:);
conl = conl(pkeep);                     % do we still need the lengths?

p = setdiff(p,shp.Points,'rows');
pfix = unique(shp.Points,'rows');
p = [pfix; p];

t=delaunayn(p);                                  % List of triangles
pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;    % Compute centroids
plot(pmid(:,1),pmid(:,2),'.r')

facereject = ~inShape(shp,pmid);
t = t(facereject,:);                                % delete interior faces
bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
bars=unique(sort(bars,2),'rows');                % Bars as node pairs

patch('vertices',p,'faces',t,'edgecol','k','facecol',[.8,.9,1]);

hold off