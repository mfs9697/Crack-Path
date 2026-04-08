function [Bound, G] = build_domain_pencil_polyline(Pmid, A, B, w, varargin)
%BUILD_DOMAIN_PENCIL_POLYLINE  Polyline pencil-channel boundary builder
%that matches the proven topology of build_polygon_pencil:
%   Bound = [ upper_chain; lower_chain; 0,-B; A,-B; A,0; A,B; 0,B ];
%and closes implicitly back to the first point.
%
% Inputs:
%   Pmid : (N x 2) midline vertices, ordered mouth -> ... -> tip
%   A,B  : rectangle size, x in [0,A], y in [-B,B]
%   w    : channel half-width
%
% Name-value options:
%   'join'        : 'miter' (default) or 'bevel'
%   'miter_limit' : positive scalar (default 6)
%   'corner_tol'  : positive scalar (default 1e-10)
%   'tip'         : 'point' (default) or 'offset'
%
% Output:
%   Bound : (M x 2) CCW polygon vertices
%   G     : debug struct with chains and landmarks

    % ---- Parse options (MATLAB 2023a friendly) ----
    p = inputParser;
    p.addParameter('join', 'miter', @(s)ischar(s)||isstring(s));
    p.addParameter('miter_limit', 6, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    p.addParameter('corner_tol', 1e-10, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    p.addParameter('tip', 'point', @(s)ischar(s)||isstring(s));
    p.parse(varargin{:});
    opts = p.Results;

    if size(Pmid,1) < 2
        error('Pmid must have at least 2 points.');
    end

    % ---- Build offset chains around polyline ----
    [Up, Dn, dbg] = offset_polyline(Pmid, w, opts);

    up_mouth = Up(1,:);
    dn_mouth = Dn(1,:);
    tip      = Pmid(end,:);

    % Snap mouth x to the left boundary (helps robustness)
    up_mouth(1) = 0;  dn_mouth(1) = 0;
    Up(1,:) = up_mouth; Dn(1,:) = dn_mouth;

    % ---- Channel boundary: upper to tip, then lower back to mouth ----
    upper_chain = Up;            % mouth -> tip

    lower_chain = flipud(Dn);    % tip -> mouth
    % Avoid duplicating the tip point if tip is apex (Up(end)=Dn(end)=tip)
    if norm(upper_chain(end,:) - lower_chain(1,:)) < 1e-14
        lower_chain = lower_chain(2:end,:);
    end

    % In your proven builder, lower_chain starts from interior points and ends at dn_mouth.
    % Here, lower_chain already ends at dn_mouth, which is what we want.

    % ---- Assemble boundary EXACTLY in the old robust topology ----
    % No explicit insertion of up_mouth at the end; closure is implicit.
    Bound = [ upper_chain;
              lower_chain;
              0, -B;
              A, -B;
              A,  0;
              A,  B;
              0,  B ];

    % Remove consecutive duplicates (prevents zero-length edges)
    Bound = Bound([true; any(diff(Bound,1,1)~=0,2)], :);

    % Ensure CCW
    if signed_area(Bound) < 0
        Bound = flipud(Bound);
    end

    % ---- Debug outputs ----
    G.up_mouth = up_mouth;
    G.dn_mouth = dn_mouth;
    G.tip      = tip;
    G.up_chain = Up;
    G.dn_chain = Dn;
    G.debug    = dbg;
end

% ========================= helpers =========================

function [Up, Dn, dbg] = offset_polyline(P, w, opts)
%OFFSET_POLYLINE  Build upper/lower offsets of a polyline with robust corner handling.

    N = size(P,1);

    % Segment tangents and normals
    T = zeros(N-1,2);
    Nrm = zeros(N-1,2);
    for i = 1:N-1
        v = P(i+1,:) - P(i,:);
        t = unit2(v);
        T(i,:) = t;
        Nrm(i,:) = [-t(2), t(1)];
    end

    Up = zeros(N,2);
    Dn = zeros(N,2);

    % Mouth endpoint: offset by first normal
    Up(1,:) = P(1,:) + w*Nrm(1,:);
    Dn(1,:) = P(1,:) - w*Nrm(1,:);

    % Tip endpoint behavior
    if strcmp(opts.tip,'point')
        Up(N,:) = P(N,:);
        Dn(N,:) = P(N,:);
    else
        Up(N,:) = P(N,:) + w*Nrm(end,:);
        Dn(N,:) = P(N,:) - w*Nrm(end,:);
    end

    % Interior vertices
    for k = 2:N-1
        t1 = T(k-1,:);
        t2 = T(k,:);
        turn = abs(t1(1)*t2(2) - t1(2)*t2(1)); % |cross2(t1,t2)|

        % Smooth continuation points
        Au = P(k,:) + w*Nrm(k,:);
        Ad = P(k,:) - w*Nrm(k,:);

        if turn < opts.corner_tol
            Up(k,:) = Au;
            Dn(k,:) = Ad;
            continue
        end

        % Corner points
        A1 = P(k,:) + w*Nrm(k-1,:);  v1 = t1;
        A2 = P(k,:) + w*Nrm(k,:);    v2 = t2;

        B1 = P(k,:) - w*Nrm(k-1,:);  u1 = t1;
        B2 = P(k,:) - w*Nrm(k,:);    u2 = t2;

        if strcmp(opts.join,'bevel')
            Up(k,:) = A2;
            Dn(k,:) = B2;
        else
            [pu, okU] = line_intersect_safe(A1, v1, A2, v2);
            [pd, okD] = line_intersect_safe(B1, u1, B2, u2);

            if ~okU || norm(pu - Au) > opts.miter_limit*w
                pu = A2;
            end
            if ~okD || norm(pd - Ad) > opts.miter_limit*w
                pd = B2;
            end

            Up(k,:) = pu;
            Dn(k,:) = pd;
        end
    end

    dbg.T = T;
    dbg.Nrm = Nrm;
end

function [p, ok] = line_intersect_safe(pA, vA, pB, vB)
    M = [vA(:), -vB(:)];
    rhs = (pB - pA).';
    if rcond(M) < 1e-12
        p = [NaN, NaN];
        ok = false;
        return;
    end
    ts = M \ rhs;
    p = pA + ts(1)*vA;
    ok = all(isfinite(p));
end

function t = unit2(v)
    nv = hypot(v(1), v(2));
    if nv < 1e-30
        t = [1,0];
    else
        t = v / nv;
    end
end

function A = signed_area(P)
    x = P(:,1); y = P(:,2);
    A = 0.5*sum(x.*y([2:end 1]) - x([2:end 1]).*y);
end
