% [TR,shp] = execAlphaShp(xyz_data(:,1:2),'Holethreshold',1e4);
[fbConn, fbFacets] = TR.freeBoundary;
bboxVert = [bbox;bbox(:,1),flipud(bbox(:,2))]; 
bboxVertConn = [1 3;3 2;2 4;4 1];

node = [fbFacets(:,1:2);bboxVert];
edge = [fbConn;bboxVertConn+length(fbConn)];


% aus tridemo von mesh2d demo 3:
%---------------------------------------------- do size-fun.
   [vlfs,tlfs, hlfs] = lfshfn2(node,edge) ;

   [slfs] = idxtri2(vlfs,tlfs) ;

%---------------------------------------------- do mesh-gen.
    hfun = @trihfn2;

    opts = struct();
    opts.debug = true;
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],opts,hfun , ...
                    vlfs,tlfs,slfs,hlfs);

%---------------------------------------------- do mesh-opt.
   [vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;

    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;

    tricost(vert,etri,tria,tnum);

    drawnow;

    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.05,.05,.30,.35]) ;