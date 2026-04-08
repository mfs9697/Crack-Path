function [mesh, mat, quad, elod, Stif, cohes] = geom_pencil(C)
%GEOM_PENCIL  Polyline-capable geometry + PDE mesh generator with pencil channel.
%
%   [V, mesh, mat, quad, elod, Stif, cohes] = geom_pencil(C)
%
% Optimized version:
%   - T3 -> T6 via T3toT6_fast (global edge map)
%   - cohesive midnodes via edge2mid map (no O(ne) findMidNode loop)
%   - last cohesive leg is straight => node picking by ONE segment per face
%   - fix t-mask ordering: mask first, then clamp
%
% Notes:
%   - returns cohes as numeric (2*K+1)x2 index list (upper/lower), interleaved v,m,v,...
%   - cohesive metadata stored in mesh.coh

        % -------------------- defaults / required --------------------
    must(C,'A'); must(C,'B'); must(C,'chw'); must(C,'ncoh');

    A    = C.A;
    B    = C.B;
    w    = C.chw;
    ncoh = max(1, round(C.ncoh));

    % Crack polyline input
    if isfield(C,'Pmid') && ~isempty(C.Pmid)
        Pmid0 = C.Pmid;
    else
        error('C must contain Pmid.');
    end

    if size(Pmid0,1) < 2
        error('C.Pmid must have at least 2 points.');
    end

    % Options
    join        = getf(C,'join','miter');      % 'miter' or 'bevel'
    miter_limit = getf(C,'miter_limit',6);
    corner_tol  = getf(C,'corner_tol',1e-10);
    tip_mode    = getf(C,'tip','point');       % 'point' recommended
    hgrad       = getf(C,'hgrad',1.15);
    hmax_ratio  = getf(C,'hmax_ratio',40);

    % Optional holes
    if ~isfield(C,'holes') || isempty(C.holes)
        C.holes = {};
    end

    % -------------------- subdivide last leg --------------------
    Pmid = subdivide_last_leg(Pmid0, ncoh);

    Llast = norm(Pmid0(end,:) - Pmid0(end-1,:));
    hcoh  = Llast / ncoh;

    % -------------------- build outer boundary polygon --------------------
    [BoundOuter, G] = build_domain_pencil_polyline_topology( ...
        Pmid, A, B, w, ...
        'join', join, ...
        'miter_limit', miter_limit, ...
        'corner_tol', corner_tol, ...
        'tip', tip_mode );

    % -------------------- build hole loops --------------------
    holeLoops = holes_to_loops(C.holes);

    % -------------------- PDE geometry + mesh --------------------
    Hcap  = min(B/2, B / hmax_ratio);
    Hmax  = max(Hcap, 1.05*hcoh);
    Hmin  = max(0.98*hcoh, 1e-12);
    Hgrad = max(1.01, hgrad);

    mdl = build_pde_geometry_with_holes(BoundOuter, holeLoops);

    msh = generateMesh(mdl, ...
        'Hmin', Hmin, ...
        'Hmax', Hmax, ...
        'Hgrad', Hgrad, ...
        'GeometricOrder','linear');
    coord0   = msh.Nodes.';     % nnode x 2
    connect0 = msh.Elements.';  % nelem x 3

    % -------------------- outputs (match “existing fields” style) --------------------
    mesh = struct();
    mesh.model    = mdl;
    mesh.meshobj  = msh;
    mesh.coord0   = coord0;
    mesh.connect0 = connect0;

    mesh.Bound = BoundOuter;
    mesh.Ggeom = G;
    mesh.Pmid0 = Pmid0;
    mesh.Pmid  = Pmid;

        % optional plots
    if getf(C,'plotGeom',true)
        figure(30); clf; hold on; axis equal; box on

        % Outer boundary
        plot([BoundOuter(:,1); BoundOuter(1,1)], ...
             [BoundOuter(:,2); BoundOuter(1,2)], ...
             'k-', 'LineWidth', 1.2);

        % Hole boundaries
        for ih = 1:numel(holeLoops)
            H = holeLoops{ih};
            plot([H(:,1); H(1,1)], ...
                 [H(:,2); H(1,2)], ...
                 'k-', 'LineWidth', 1.2);
        end

        % Crack midline
        plot(Pmid(:,1), Pmid(:,2), 'r.-', 'LineWidth', 1.4, 'MarkerSize', 12);

        % Offset chains
        plot(G.up_chain(:,1), G.up_chain(:,2), 'b--');
        plot(G.dn_chain(:,1), G.dn_chain(:,2), 'b--');

        title('Geometry: outer boundary + hole(s) + midline + offsets');
        xlim([0 A]);
        ylim([-B B]);
    end

    if getf(C,'plotMesh',true)
        figure(31); clf; hold on; axis equal; box on

        % Linear mesh
        triplot(connect0, coord0(:,1), coord0(:,2));

        % Outer boundary
        plot([BoundOuter(:,1); BoundOuter(1,1)], ...
             [BoundOuter(:,2); BoundOuter(1,2)], ...
             'r-', 'LineWidth', 1.0);

        % Hole boundaries
        for ih = 1:numel(holeLoops)
            H = holeLoops{ih};
            plot([H(:,1); H(1,1)], ...
                 [H(:,2); H(1,2)], ...
                 'r-', 'LineWidth', 1.0);
        end

        % Crack midline
        plot(Pmid(:,1), Pmid(:,2), 'k.-', 'LineWidth', 1.0, 'MarkerSize', 10);

        title('Mesh: T3 + outer boundary + hole(s) + midline');
        xlim([0 A]);
        ylim([-B B]);
    end
    % ========================= Collapse physical crack + cohesive leg =========================
    coord_lin   = coord0;
    connect_lin = connect0;

    leave_tip_gap = true;
    if isfield(C,'leave_tip_gap'), leave_tip_gap = logical(C.leave_tip_gap); end

    % Original midline polyline (mouth -> ... -> knee -> tip)
    P0 = mesh.Pmid0;
    Np = size(P0,1);
    if Np < 2
        error('mesh.Pmid0 must contain at least 2 points.');
    end

    % Cohesive leg is the last segment (knee -> tip)
    p1 = P0(end-1,:);          % knee
    p2 = P0(end,:);            % tip

    % Where cohesive leg begins in subdivided polyline
    Psub = mesh.Pmid;
    Nsub = size(Psub,1);
    k0   = Nsub - ncoh;        % index of p1 in mesh.Pmid

    % Channel face chains (built from offset_polyline on mesh.Pmid)
    Up_chain = mesh.Ggeom.up_chain;
    Dn_chain = mesh.Ggeom.dn_chain;

    % Tolerance for selecting nodes near channel faces (works even if only endpoints are boundary vertices)
    tol_line = max(1e-3*Hmin, 1e-6*w);

    % -------------------------------------------------------------------------
    % (A) Physical crack collapse (all segments except the last): projection to midline
    % -------------------------------------------------------------------------
    % Physical part in P0 is segments 1..(Np-2) (since last segment is cohesive)
    if Np > 2
        for i = 1:(Np-2)
            % midline segment
            a = P0(i,:);
            b = P0(i+1,:);
            AB = b - a;
            L  = norm(AB);
            if L < 1e-30
                continue
            end
            e  = AB / L;

            % corresponding channel-face segments (same index i->i+1)
            Au = Up_chain(i,:);   Bu = Up_chain(i+1,:);
            Ad = Dn_chain(i,:);   Bd = Dn_chain(i+1,:);

            % select ALL mesh nodes close to the upper/lower face segments
            iu = face_ids_on_segment(coord0, Au, Bu, tol_line);
            id = face_ids_on_segment(coord0, Ad, Bd, tol_line);

            ids = unique([iu; id], 'stable');
            if isempty(ids)
                continue
            end

            % project to the midline segment a->b (no midpoint averaging!)
            t = (coord_lin(ids,:) - a) * e.';   % arclength along segment
            t = max(0, min(L, t));
            coord_lin(ids,:) = a + t.*e;
        end
    end

    % -------------------------------------------------------------------------
    % (B) Cohesive leg collapse + indices (keep your existing approach)
    % -------------------------------------------------------------------------
    % cohesive faces on channel boundary (knee -> ... -> tip) from subdivided chains
    Up2 = Up_chain(k0:end,:);
    Dn2 = Dn_chain(k0:end,:);

    [leg2_up, sU] = face_ids_on_polyline_param(coord0, Up2, tol_line);
    [leg2_dn, sD] = face_ids_on_polyline_param(coord0, Dn2, tol_line);

    if isempty(leg2_up) || isempty(leg2_dn)
        error('Cohesive face node sets are empty. Increase chw or relax tol_line.');
    end

    ttar = linspace(0,1,ncoh+1).';
    idsU = pick_unique_by_targets(sU, leg2_up, ttar);
    idsD = pick_unique_by_targets(sD, leg2_dn, ttar);
    pairs_v = [idsU(:), idsD(:)];     % (ncoh+1) x 2

    % collapse to midpoints (optionally keep last pair as tip gap)
    if leave_tip_gap && size(pairs_v,1) > 1
        idx = 1:(size(pairs_v,1)-1);
    else
        idx = 1:size(pairs_v,1);
    end

    Cmid = 0.5*(coord_lin(pairs_v(idx,1),:) + coord_lin(pairs_v(idx,2),:));
    coord_lin(pairs_v(idx,1),:) = Cmid;
    coord_lin(pairs_v(idx,2),:) = Cmid;

    % snap knee pair exactly to p1
    coord_lin(pairs_v(1,1),:) = p1;
    coord_lin(pairs_v(1,2),:) = p1;

    % clean degenerate triangles (after collapse) and enforce CCW
    Atri = triAreasSigned(connect_lin, coord_lin);
    keep = abs(Atri) > (1e-12 * max(hcoh,1e-12)^2);  % scaled threshold
    connect_lin = connect_lin(keep,:);
    Atri = Atri(keep);

    cw = Atri < 0;
    if any(cw)
        tmp = connect_lin(cw,2);
        connect_lin(cw,2) = connect_lin(cw,3);
        connect_lin(cw,3) = tmp;
    end

    % ---------- OPTIMIZATION 3: upgrade to T6 using fast global edge map ----------
    [coord, connect] = T3toT6_fast(coord_lin, connect_lin);

    %{
    figure; hold on; axis equal; box on

    nt = 12;                         % resolution per edge
    t  = linspace(0,1,nt).';

    for e = 1:size(connect,1)
        v = connect(e,1:3);
        m = connect(e,4:6);

        % edges: (v1,m12,v2), (v2,m23,v3), (v3,m31,v1)
        edges = [
            v(1), m(1), v(2);
            v(2), m(2), v(3);
            v(3), m(3), v(1)
            ];

        Px = []; Py = [];

        for k = 1:3
            P = coord(edges(k,:),:);   % [v_i; m_ij; v_j]

            B = (1-t).^2 .* P(1,:) ...
                + 2*(1-t).*t .* P(2,:) ...
                + t.^2 .* P(3,:);

            Px = [Px; B(:,1)];
            Py = [Py; B(:,2)];
        end

        patch(Px, Py, [0.9 0.9 0.9], ...
            'EdgeColor','k', ...
            'LineWidth',0.8);
    end

    title('T6 mesh: true quadratic element patches');
    xlabel('x'); ylabel('y');

    %}

    % ---------- OPTIMIZATION 4: build mid-side cohesive pairs using edge2mid map ----------
    K = size(pairs_v,1) - 1;
    if K < 1
        cohes = zeros(0,2);
        pairs_m = zeros(0,2);
    else
        edge2mid = build_edge2mid_map_T6(connect);

        pairs_m = zeros(K,2);
        for i = 1:K
            pairs_m(i,1) = lookup_midnode(edge2mid, pairs_v(i,1), pairs_v(i+1,1));
            pairs_m(i,2) = lookup_midnode(edge2mid, pairs_v(i,2), pairs_v(i+1,2));
        end

        % Interleave: v1, m1, v2, m2, ..., v_{K}, m_{K}, v_{K+1}
        cohes = zeros(2*K+1, 2);
        cohes(1:2:end,:) = pairs_v;
        cohes(2:2:end,:) = pairs_m;
    end

    % store final mesh
    mesh.coord    = coord;
    mesh.connect  = connect;

    % cohesive bookkeeping
    mesh.coh = struct();
    mesh.coh.p1 = p1; mesh.coh.p2 = p2;
    mesh.coh.targets = ttar;
    mesh.coh.pairs_v = pairs_v;
    mesh.coh.pairs_m = pairs_m;

    % -------------------- physics / loads / stiffness --------------------
    elod = edge_loads_T6(coord, B, C.eps1);

    mat  = struct('E',C.E2, 'nu',C.nu, 'G12',C.G12, 'D',C.Dmat);

    [nip2, xip2, w2, Nextr] = integr();
    quad = struct('nip2',nip2,'xip2',xip2,'w2',w2,'Nextr',Nextr);

    % essential BCs (remove rigid modes)
    fix_pts = [A,0; 0,B; 0,-B];
    fix = zeros(3,1);
    for i=1:3
        [~,fix(i)] = min((coord(:,1)-fix_pts(i,1)).^2 + (coord(:,2)-fix_pts(i,2)).^2);
    end
    if numel(unique(fix)) < 3
        warning('BC points collided (mesh too coarse). Consider refining Hmax or adjusting fix_pts.');
    end
    fixvar = [2*fix(1); 2*fix(2)-1; 2*fix(3)-1];

    Stif = stif_assem(mesh, mat, quad, fixvar);

    % ========================= helper functions =========================
    function must(C, field)
        if ~isfield(C,field)
            error('Missing required field C.%s', field);
        end
    end

    function v = getf(C, field, default)
        if isfield(C,field) && ~isempty(C.(field))
            v = C.(field);
        else
            v = default;
        end
    end

    function ids = pick_unique_by_targets(tall, idsAll, ttar)
        % tall sorted ascending, idsAll aligned with tall
        % returns one unique id per target, monotone.
        n = numel(tall);
        m = numel(ttar);
        ids = zeros(m,1);

        i = 1;
        for k = 1:m
            tk = ttar(k);
            if i > n
                error('pick_unique_by_targets: not enough face nodes to satisfy ncoh+1 targets.');
            end

            best_i = i;
            best_d = abs(tall(i) - tk);

            while (best_i + 1) <= n
                dnext = abs(tall(best_i+1) - tk);
                if dnext <= best_d
                    best_i = best_i + 1;
                    best_d = dnext;
                else
                    break;
                end
            end

            ids(k) = idsAll(best_i);
            i = best_i + 1; % enforce uniqueness
        end
    end

    function ids = face_ids_on_segment(P, A, B, tol)
        if nargin < 4 || isempty(tol), tol = 1e-12; end
        AB = B - A; L2 = max(1e-32, dot(AB,AB));
        n = [-AB(2), AB(1)];
        nn = norm(n);
        if nn==0
            d = zeros(size(P,1),1);
        else
            d = (P-A)*(n.'/nn);
        end
        t = ((P - A) * AB.') / L2;
        mask = (abs(d) < tol) & (t >= -1e-10) & (t <= 1+1e-10);
        ids  = find(mask);
        if isempty(ids), return, end
        s = vecnorm(P(ids,:) - A, 2, 2);
        ids = sortrows([ids, s], 2); ids = ids(:,1);
    end

    function [ids, s] = face_ids_on_polyline_param(P, Q, tol)
        %FACE_IDS_ON_POLYLINE_PARAM  Collect node ids near polyline Q and return
        %arc-length parameter s in [0,1] for sorting/targeting.
        %
        % P : [N x 2] node coordinates
        % Q : [M x 2] polyline points
        % tol : distance tolerance
        %
        % ids : indices into P
        % s   : corresponding arc-length fraction along Q (0 at Q(1), 1 at Q(end))

        M = size(Q,1);
        if M < 2
            ids = []; s = [];
            return;
        end

        % bounding-box prefilter
        xmin = min(Q(:,1)) - tol; xmax = max(Q(:,1)) + tol;
        ymin = min(Q(:,2)) - tol; ymax = max(Q(:,2)) + tol;
        cand = find(P(:,1) >= xmin & P(:,1) <= xmax & P(:,2) >= ymin & P(:,2) <= ymax);

        if isempty(cand)
            ids = []; s = [];
            return;
        end

        X = P(cand,:);

        % segments
        P1 = Q(1:end-1,:);
        P2 = Q(2:end,:);
        V =  P2 - P1;  % [M-1 x 2]
        L = hypot(V(:,1), V(:,2));      % segment lengths
        L(L < 1e-30) = 1e-30;
        Vn = V ./ L;                    % unit tangents

        % cumulative arc length
        cumL = [0; cumsum(L)];
        Ltot = cumL(end);

        % find closest segment to each X
        d2_best = inf(size(X,1),1);
        s_best  = zeros(size(X,1),1);

        for k = 1:(M-1)
            Ak = P1(k,:); vk = Vn(k,:); lk = L(k);

            % projection parameter on segment [0, lk]
            r  = X - Ak;                        % [Nc x 2]
            tk = r*vk.';                        % [Nc x 1] signed distance along segment
            tk = max(0, min(lk, tk));           % clamp to segment

            % closest point and distance^2
            Ck = Ak + tk*vk;                    % [Nc x 2]
            dk = X - Ck;
            d2 = dk(:,1).^2 + dk(:,2).^2;

            better = d2 < d2_best;
            if any(better)
                d2_best(better) = d2(better);
                % arc-length position along whole polyline
                s_best(better) = (cumL(k) + tk(better)) / max(Ltot,1e-30);
            end
        end

        keep = d2_best <= tol^2;
        ids  = cand(keep);
        s    = s_best(keep);

        % sort by s and unique (stable)
        [s, I] = sort(s);
        ids = ids(I);

        [ids, ia] = unique(ids, 'stable');
        s = s(ia);
    end


    function A = triAreasSigned(T, X)
        v1=X(T(:,1),:); v2=X(T(:,2),:); v3=X(T(:,3),:);
        A  = 0.5*((v2(:,1)-v1(:,1)).*(v3(:,2)-v1(:,2)) - (v2(:,2)-v1(:,2)).*(v3(:,1)-v1(:,1)));
    end

    function [vert6, tria6] = T3toT6_fast(vert, tria)
        % Upgrade T3 -> T6 using a global edge map
        ne = size(tria,1);
        n0 = size(vert,1);

        E12 = sort(tria(:,[1 2]),2);
        E23 = sort(tria(:,[2 3]),2);
        E31 = sort(tria(:,[3 1]),2);

        Eall = [E12; E23; E31];
        [Eu, ~, ic] = unique(Eall, 'rows');
        nu = size(Eu,1);

        mids = (n0+1 : n0+nu).';
        vert6 = [vert; 0.5*(vert(Eu(:,1),:) + vert(Eu(:,2),:))];

        m12 = mids(ic(1:ne));
        m23 = mids(ic(ne+1:2*ne));
        m31 = mids(ic(2*ne+1:3*ne));

        tria6 = [tria, m12, m23, m31];
    end

    function edge2mid = build_edge2mid_map_T6(conn6)
        v = conn6(:,1:3);
        m = conn6(:,4:6);

        E12 = sort([v(:,1), v(:,2)],2);  M12 = m(:,1);
        E23 = sort([v(:,2), v(:,3)],2);  M23 = m(:,2);
        E31 = sort([v(:,3), v(:,1)],2);  M31 = m(:,3);

        E = [E12; E23; E31];
        M = [M12; M23; M31];

        [Eu, ia] = unique(E, 'rows', 'stable');
        Mu = M(ia);

        edge2mid.E = Eu;
        edge2mid.M = Mu;
    end

    function mid = lookup_midnode(edge2mid, v1, v2)
        e = sort([v1 v2]);
        [tf, loc] = ismember(e, edge2mid.E, 'rows');
        if ~tf
            error('Mid node for edge (%d,%d) not found.', v1, v2);
        end
        mid = edge2mid.M(loc);
    end
end

