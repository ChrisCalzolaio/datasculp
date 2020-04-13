function [p,t]=distmesh2dDiscrete(shp,h0,bbox,figH,varargin)
%DISTMESH2D 2-D Mesh Generator using Distance Functions.
%   [P,T]=DISTMESH2D(FD,FH,H0,BBOX,PFIX,FPARAMS)
%{
%      P:         Node positions (Nx2)
%      T:         Triangle indices (NTx3)
%      FD:        Distance function d(x,y)
%      FH:        Scaled edge length function h(x,y)
%      H0:        Initial edge length
%      BBOX:      Bounding box [xmin,ymin; xmax,ymax]
%      PFIX:      Fixed node positions (NFIXx2)
%      FPARAMS:   Additional parameters passed to FD and FH
%
%       from: https://people.sc.fsu.edu/~jburkardt/m_src/distmesh/distmesh.html
%           fd, the name of a distance function defining the region;
%           fh, the name of a mesh density function;
%           h, the nominal mesh spacing;
%           box, defining a box that contains the region;
%           iteration_max, the maximum number of iterations;
%           fixed, a list of points which must be included in the mesh, or '[]', if no fixed points are given.
%           p, (output), a list of node coordinates;
%           t, (output), a list of node indices forming triangles;
%}
%{
%   See also: MESHDEMO2D, DISTMESHND, DELAUNAYN, TRIMESH.
%   distmesh2d [Discrete].m v1.1
%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
%}

% graphical output:
for fig=1:numel(figH)
    set(0,'CurrentFigure',figH(fig));
    clf;
end

% logging
varnames = {'Step','Duration','DeltaPointsAV','NumberPoints','st3t','st4t','st5t','st6t'};
vartypes = repmat({'double'},1,numel(varnames));
logt = table('Size',[0 numel(varnames)],'VariableTypes',vartypes,'VariableNames',varnames);
clearvars varnames vartypes;
% config
dptol   = 0.001;        % delta points tolerance, used for termination criterion, original value: 0.001
ttol    = 0.1;
Fscale  = 1.2;
deltat  = 0.2;          % original value; 0.2
geps    = 0.001*h0;
deps    = sqrt(eps)*h0;
densityctrlfreq=30;

bboxoffset = diff(bbox)/2;
bboxoffset = [-bboxoffset;bboxoffset];
extbbox = bbox + bboxoffset;  % extended bounding box

set(0,'CurrentFigure',figH(2));
plot(shp,'LineStyle','none');

% 1. Create initial distribution in bounding box (equilateral triangles)
tic
[x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
x(2:2:end,:)=x(2:2:end,:)+h0/2;                      % Shift even rows
p=[x(:),y(:)];                                       % List of node coordinates
line(p(:,1),p(:,2),'Marker','.','Color','y','LineStyle','none');
fprintf('[ %s ] distmesh2d [Discrete]: step 1 [ %.3f ]\n',datestr(now,'HH:mm:ss'),toc)

% 2. Remove points outside the region, apply the rejection method
tic;
p = p(~shp.inShape(p),:);                                                   % keep points in p outside of shape
% line(p(:,1),p(:,2),'Marker','.','Color','b','LineStyle','none');
[~,distBoundary] = shp.nearestNeighbor(p);                                  % get absolute distance from boundary
p = p(rand(size(p,1),1)>distBoundary ./ max(distBoundary),:);               % keep only points with rand number greater than normalized distance from boundary

[~,pBoundary] = shp.boundaryFacets;                                         % find points on the boundary of the alpha shape
bboxVert = [bbox;bbox(:,1),flipud(bbox(:,2))];                              % fixed points at the corners of the bounding box
% pFixed = [pBoundary; bboxVert];
pFixed = pBoundary;
nfix=length(pFixed);                                                        % number of fixed points

p = setdiff(p,pFixed,'rows');

line(p(:,1),p(:,2),'Marker','.','Color','b','LineStyle','none');

pfix = logical([zeros(length(p),1);ones(nfix,1)]);                          % keep metadata which points are fixed (fixed=true)
p = [p; pFixed];                                               % merge outside points with boundary points
N=length(p);                                                                % Number of points N

% p=p(feval(fd,p,varargin{:})<geps,:);                 % Keep only d<0 points, keeps points within the shape
% r0=1./feval(fh,p,varargin{:}).^2;                    % Probability to keep point
% p=p(rand(size(p,1),1)<r0./max(r0),:);                % Rejection method
% if ~isempty(pfix)
%     p=setdiff(p,pfix,'rows');                        % Remove duplicated nodes
% end
% pfix=unique(pfix,'rows');

% p=[pfix; p];                                         % Prepend fix points

count=0;
pold=inf;                                            % For first iteration
clf
view(2)
axis equal
% axis off
fprintf('[ %s ] distmesh2d [Discrete]: step 2 [ %.3f ]\n',datestr(now,'HH:mm:ss'),toc)
while 1
    ittime = tic;
    count=count+1;
    %     logt.iteration =
    fprintf('[ %s ] distmesh2d [Discrete]: iteration %d started.\n',datestr(now,'HH:mm:ss'),count)
    % 3. Retriangulation by the Delaunay algorithm
    st3t = tic;
    if max(sqrt(sum((p-pold).^2,2))/h0)>ttol           % Any large movement?
        set(0,'CurrentFigure',figH(2));
        pold=p;                                          % Save current positions
        t=delaunay(p);                                  % List of triangles
        triplot(t,p(:,1),p(:,2))
        pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;    % Compute centroids
        line(pmid(:,1),pmid(:,2),'Marker','.','Color','r','LineStyle','none');
        t = t(~shp.inShape(pmid),:);                    % keep exterior triangles
        triplot(t,p(:,1),p(:,2))
        %         t=t(feval(fd,pmid,varargin{:})<-geps,:);         % Keep interior triangles
        fprintf('[ %s ] distmesh2d [Discrete]: step 3 [ %.3f ]\n',datestr(now,'HH:mm:ss'),toc(st3t))
        % 4. Describe each bar by a unique pair of nodes
        st4t = tic;
        bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
        bars=unique(sort(bars,2),'rows');                % Bars as node pairs
        fprintf('[ %s ] distmesh2d [Discrete]: step 4 [ %.3f ]\n',datestr(now,'HH:mm:ss'),toc(st4t))
        % 5. Graphical output of the current mesh
        st5t = tic;
        %         cla
        %         patch('vertices',p,'faces',t,'edgecol','k','facecol',[.8,.9,1]);
        for figs=2:3
            set(0,'CurrentFigure',figH(figs));
            triplot(t,p(:,1),p(:,2))
            switch figs
                case 2
                    xline(bbox(1,1));xline(bbox(2,1));
                    yline(bbox(1,2));yline(bbox(2,2));
                    set(gca,'XLim',[-150 350],'DataAspectRatio',[1 1 1])
                case  3
                    set(gca,'XLim',[-25 60], 'YLim',[1980 2040])
            end
            drawnow
        end        
        fprintf('[ %s ] distmesh2d [Discrete]: step 5 [ %.3f ]\n',datestr(now,'HH:mm:ss'),toc(st5t))
    end
    if count == 1
        %         TR = triangulation(t,p);
        %         [~,bboxPoints] = TR.freeBoundary;
        %         bboxPoints = setdiff(bboxPoints,pBoundary,'rows');
        %         line(bboxPoints(:,1),bboxPoints(:,2),'Marker','.','Color','g','LineStyle','none');
        % find corner points
        %         NEcornerID = TR.nearestNeighbor(1100,2400)
    end
    
    % 6. Move mesh points based on bar lengths L and forces F
    st6t = tic;
    barvec=p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
    barsmid =  (p(bars(:,1),:) + p(bars(:,2),:))/2;
%     barmidscatterH = line(barsmid(:,1),barsmid(:,2),'Marker','.','Color','m','LineStyle','none');
    
    L=sqrt(sum(barvec.^2,2));                          % L = Bar lengths
    [~,distBarsMid] = shp.nearestNeighbor(barsmid);
    hbars = distBarsMid./max(distBarsMid);              % weighting of bar midpoint to boundary;
    
%     figHmax = figure();
%     figHmax.Name = 'max';
%     figHmin = figure();
%     figHmin.Name = 'min';
%     
%     set(0,'CurrentFigure',figHmax)
%     figHmaxAx = axes();
%     pcshow([barsmid,hbars])
%     set(gca,'DataAspectRatio',[1 1 100])
%     
%     set(0,'CurrentFigure',figHmin)
%     figHminAx = axes();
%     hbars = distBarsMid./min(distBarsMid);
%     pcshow([barsmid,hbars])
%     
    
    L0=hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));     % L0 = Desired lengths
    
%     Density control - remove points that are too close
        if mod(count,densityctrlfreq)==0 & any(L0>2*L)
            p(setdiff(reshape(bars(L0>2*L,:),[],1),1:nfix),:)=[];
            N=size(p,1);
            pold=inf;
            continue;
        end
    
    F=max(L0-L,0);                                     % Bar forces (scalars)
    Fvec=F./L*[1,1].*barvec;                           % Bar forces (x,y components)
    Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
    Ftot(pfix,:) = 0;                          % Force = 0 at fixed points
    %     Ftot(xmask,1) = 0;
    %     Ftot(ymask,2) = 0;
    fprintf('[ %s ] distmesh2d [Discrete]: number of points: \t%i\n',datestr(now,'HH:mm:ss'),length(p));
    %     fprintf('[ %s ] distmesh2d [Discrete]: length force array [ %i ]\n',datestr(now,'HH:mm:ss'),length(Ftot));
    if ~(length(p)==length(Ftot))
        pause();
    end
    p=p+deltat*Ftot;                                   % Update node positions
    
    fprintf('[ %s ] distmesh2d [Discrete]: step 6 [ %.3f ]\n',datestr(now,'HH:mm:ss'),toc(st6t))
    
    % 7. Bring outside points back to the boundary
    tic
%     xmask = ~(and(p(:,1)>bbox(1,1),p(:,1)<bbox(2,1)));
%     ymask = ~(and(p(:,2)>bbox(1,2),p(:,2)<bbox(2,2)));
    xmasklo = p(:,1)<bbox(1,1);
    xmaskup = p(:,1)>bbox(2,1);
    ymasklo = p(:,2)<bbox(1,2);
    ymaskup = p(:,2)>bbox(2,2);
    
    p(xmasklo,1) = bbox(1,1);
    p(xmaskup,1) = bbox(2,1);
    p(ymasklo,2) = bbox(1,2);
    p(ymaskup,2) = bbox(2,2);
%     
%     d=feval(fd,p,varargin{:});                         % Find points outside (d>0)
%     ix=d>0;
%     dgradx=(feval(fd,[p(ix,1)+deps,p(ix,2)],varargin{:})-d(ix))/deps; % Numerical
%     dgrady=(feval(fd,[p(ix,1),p(ix,2)+deps],varargin{:})-d(ix))/deps; %    gradient
%     dgrad2=dgradx.^2+dgrady.^2;
%     p(ix,:)=p(ix,:)-[d(ix).*dgradx./dgrad2,d(ix).*dgrady./dgrad2];    % Project
%     fprintf('[ %s ] distmesh2d [Discrete]: step 7 [ %.3f ]\n',datestr(now,'HH:mm:ss'),toc)
%     
    
    % 8. Termination criterion: All interior nodes move less than dptol (scaled)
    %     dp_av = max(sqrt(sum(deltat*Ftot(d<-geps,:).^2,2))/h0);             % delta point, actual value
    dp_av = max(sqrt(sum(deltat*Ftot.^2,2))/h0);
    fprintf('[ %s ] distmesh2d [Discrete]: delta points, actual value: %.5f\n',datestr(now,'HH:mm:ss'),dp_av)
    if dp_av<dptol
        break;
    end
    
    N=length(p);
    
    fprintf('[ %s ] distmesh2d [Discrete]: iteration %d done: [ %.3f ]\n',datestr(now,'HH:mm:ss'),count, toc(ittime))
    logt(count,:) = {count, ittime, dp_av, N, st3t, st4t, st5t, st6t};
    % plot
    set(0,'CurrentFigure',figH(1))
    subplot(2,1,1)
    line(logt.Step,logt.NumberPoints);
    subplot(2,1,2)
    line(logt.Step,logt.DeltaPointsAV);
    
%     plotNy(logt.Step,...
%             {logt.Duration,logt.DeltaPointsAV,logt.NumberPoints,...
%             logt.st3t,logt.st4t,logt.st5t,logt.st6t},...
%             [1 2 3 4 4 4 4],...
%             'YAxisLabels',{'time [s]','distance','#','time [s]'},...
%             'TitleStr','logging',...
%             'LegendLoc','none',...
%             'Parent',figH(3));
        drawnow;
end

% Clean up and plot final mesh
% [p,t]=fixmesh(p,t);
% simpplot(p,t)