function out = mesh_pencil_domain(C)
%MESH_PENCIL_DOMAIN  Geometry + mesh only (no physics).
%
% Uses canonical crack polyline C.Pmid (mouth -> ... -> tip), subdivides the
% last leg into C.ncoh segments, builds a pencil-like channel embedded in a
% rectangle [0,A]x[-B,B], and generates a T3 mesh (PDE toolbox).
%
% Required fields in C:
%   A, B, chw, ncoh, Pmid
% Optional fields:
%   join, miter_limit, corner_tol, tip, hmax, hgrad, plotMesh, plotGeom

    % ---- Required ----
    A   = C.A;
    B   = C.B;
    w   = C.chw;
    P0  = C.Pmid;

    if size(P0,1) < 2
        error('C.Pmid must contain at least 2 points.');
    end

    ncoh = max(1, round(C.ncoh));

    % ---- Options with defaults ----
    if ~isfield(C,'join'),        C.join = 'miter'; end
    if ~isfield(C,'miter_limit'), C.miter_limit = 6; end
    if ~isfield(C,'corner_tol'),  C.corner_tol = 1e-10; end
    if ~isfield(C,'tip'),         C.tip = 'point'; end

    if ~isfield(C,'hgrad'),       C.hgrad = 1.15; end
    if ~isfield(C,'hmax'),        C.hmax  = 0.10; end

    if ~isfield(C,'plotMesh'),    C.plotMesh = true; end
    if ~isfield(C,'plotGeom'),    C.plotGeom = true; end

    % ---- Subdivide last leg (cohesive resolution cue) ----
    Pmid = subdivide_last_leg(P0, ncoh);

    % Segment length on the (original) last leg for Hmin target
    lastSegLen = norm(P0(end,:) - P0(end-1,:));
    hcoh = lastSegLen / ncoh;

    % ---- Build domain boundary with pencil channel ----
    [Bound, Ggeom] = build_domain_pencil_polyline( ...
        Pmid, A, B, w, ...
        'join', C.join, ...
        'miter_limit', C.miter_limit, ...
        'corner_tol', C.corner_tol, ...
        'tip', C.tip );

    % ---- Visualize geometry steps (optional) ----
    if C.plotGeom
        figure(20); clf; hold on; axis equal; box on
        plot([Bound(:,1); Bound(1,1)], [Bound(:,2); Bound(1,2)], 'k-', 'LineWidth', 1.2);
        plot(Pmid(:,1), Pmid(:,2), 'r.-', 'LineWidth', 1.5, 'MarkerSize', 12);
        plot(Ggeom.up_chain(:,1), Ggeom.up_chain(:,2), 'b--');
        plot(Ggeom.dn_chain(:,1), Ggeom.dn_chain(:,2), 'b--');
        legend('Bound (domain+slot)','Midline (subdivided)','Offset +w','Offset -w', ...
            'Location','bestoutside');
        title('Geometry: domain boundary with pencil channel');
        xlim([0 A]); ylim([-B B]);
    end

    

    % ---- PDE geometry ----
    Gcol = [2; size(Bound,1); Bound(:,1); Bound(:,2)];
    gd   = decsg(Gcol);

    mdl = createpde();
    geometryFromEdges(mdl, gd);

    %{
    figure(2); clf(2); 
    pdegplot(mdl); hold on
    for i=1:size(Bound,1)
        plot(Bound(i,1),Bound(i,2),'bo'); 
    end
    %}

    % ---- Mesh ----
    delta = norm(C.Pmid(end,:)-Pmid(end-1,:));
    Hmin = max(delta/C.ncoh, 1e-12);
    Hmax = min(C.B/2, C.B/C.hmax_ratio);
    Hgrad = max(1.01, C.hgrad);

    msh = generateMesh(mdl, ...
        'Hmin', Hmin, ...
        'Hmax', Hmax, ...
        'Hgrad', Hgrad, ...
        'GeometricOrder', 'linear');

    coord   = msh.Nodes.';     % (nnode x 2)
    connect = msh.Elements.';  % (nelem x 3)

    % ---- Mesh plot (optional) ----
    if C.plotMesh
        figure(21); clf; hold on; axis equal; box on
        triplot(connect, coord(:,1), coord(:,2));
        plot([Bound(:,1); Bound(1,1)], [Bound(:,2); Bound(1,2)], 'r-', 'LineWidth', 1.1);
        plot(Pmid(:,1), Pmid(:,2), 'k.-', 'LineWidth', 1.2, 'MarkerSize', 10);
        title('Mesh + boundary + midline');
        xlim([0 A]); ylim([-B B]);
    end

    % ---- Outputs ----
    out = struct();
    out.Pmid0   = P0;        % original vertices
    out.Pmid    = Pmid;      % subdivided polyline
    out.Bound   = Bound;     % domain boundary polygon passed to decsg
    out.Ggeom   = Ggeom;     % geometry debug info (offset chains, etc.)
    out.model   = mdl;       % PDE model (geometry)
    out.mesh    = msh;       % PDE mesh object
    out.coord   = coord;     % nodes (n x 2)
    out.connect = connect;   % elements (m x 3)

    out.hcoh    = hcoh;
    out.Hmin    = Hmin;
    out.Hmax    = Hmax;
    out.Hgrad   = Hgrad;
end
