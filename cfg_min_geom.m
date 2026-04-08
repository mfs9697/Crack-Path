function C = cfg_min_geom(theta_k, varargin)
%CFG_MIN_GEOM  Minimal geometry-only configuration for building the
%domain boundary with a pencil-like channel cut.
%
% Canonical representation:
%   C.Pmid   � polyline vertices (mouth -> ... -> tip), size (nLeg+1)x2
%   C.L      � leg lengths, size nLegx1
%   C.theta  � leg global angles (rad), size nLegx1
%   C.nLeg   � number of legs
%
% Sweep usage (coarse sweep over the last leg angle):
%   C = cfg_min_geom(theta2);
%
% Extension to multi-leg trajectories:
%   Provide C.Lphys and C.thetaphys externally, then set trial cohesive leg
%   as last entry. This cfg keeps a clean, scalable representation.
%
% Optional name-value overrides:
%   'P0'      : mouth point (1x2)
%   'L'       : leg lengths (nLegx1)  (overrides default)
%   'theta'   : leg angles  (nLegx1)  (overrides default)
%   'ncoh'    : cohesive discretization count for the LAST leg
%   'hmax_ratio','hgrad','A','B','join','miter_limit','corner_tol','tip', etc.
%
% Used by:
%   [Bound,G] = build_domain_pencil_polyline(C.Pmid, C.A, C.B, C.chw);

       % ----------------- defaults -----------------
    if nargin < 1 || isempty(theta_k)
        theta_k = 0*pi/180;   % horizontal last-leg angle (rad)
    end

    % ---- Outer domain (rectangle) ----
    C.A = 0.3;     % x in [0, A]
    C.B = 0.2;     % y in [-B, B]

    % ---- Default 2-leg crack: physical + cohesive ----
    % Crack is taken horizontal along y = 0
    a      = 0.05;        % physical crack length
    da     = 0.08*a;      % cohesive probe length
    theta1 = 0*pi/180;    % first leg horizontal

    % initial window for theta
    C.init_left_deg = 1;
    C.init_right_deg = 1;

    % default bounds for adaptive delta (fictitious crack length)
    C.deltaMinFactor = 0.7;   % lower bound for trial cohesive length: delta >= deltaMinFactor * delta0
    C.deltaMaxFactor = 1.5;   % upper bound for trial cohesive length: delta <= deltaMaxFactor * delta0
    C.deltaFactor    = 1.15;   % target scaling: delta_target = deltaFactor * ell_adv (keeps active zone within trial leg)
    C.deltaRelax     = 0.5;   % relaxation parameter for update: delta_{k+1} = (1-relax)*delta_k + relax*delta_target

    C.P0    = [0.00, 0.00];   % mouth point
    C.L     = [a; da];
    C.theta = [theta1; theta_k];
    C.nLeg  = numel(C.L);

    % Current tip of the physical crack:
    x_tip = C.P0(1) + a*cos(theta1);
    y_tip = C.P0(2) + a*sin(theta1);

    % ---- Hole below and ahead of the current crack tip ----
    

    % ---- One hole ----
    %{
    % Place the hole:
    %   - ahead of the tip in +x
    %   - below the crack line in -y
    
    dx_hole = 0.06;   % horizontal shift ahead of tip
    dy_hole = 0.06;   % vertical shift downward
    R_hole  = 0.03;   % hole radius

    C.holes = {
        struct( ...
            'type',   'circle', ...
            'center', [x_tip + dx_hole, y_tip - dy_hole], ...
            'r',      R_hole, ...
            'npoly',  60 )
    };
    %}
    
    % ---- Two holes: vertically staggered ----
    % --- Hole 1: below and closer (dominant early influence) ---
    dx1 = 0.05;     % closer to tip
    dy1 = 0.05;     % downward offset
    R1  = 0.025;    % slightly smaller

    % --- Hole 2: above and farther (delayed influence) ---
    dx2 = 0.12;     % farther downstream
    dy2 = 0.05;     % upward offset
    R2  = 0.03;     % slightly larger to compensate distance

    C.holes = {
        struct( ...
            'type',   'circle', ...
            'center', [x_tip + dx1, y_tip - dy1], ...
            'r',      R1, ...
            'npoly',  60 ), ...
        struct( ...
            'type',   'circle', ...
            'center', [x_tip + dx2, y_tip + dy2], ...
            'r',      R2, ...
            'npoly',  60 )
    };

    % ---- Cohesive discretization (LAST leg) ----
    C.ncoh       = 30;
    C.hmax_ratio = 30;
    C.hgrad      = 1.1;

    % ---- Polygon join options ----
    C.join        = 'miter';     % 'miter' or 'bevel'
    C.miter_limit = 6;           % cap on miter length / chw

    % ---- Corner detection tolerance ----
    C.corner_tol  = 1e-10;

    % ---- Pencil tip behavior ----
    C.tip         = 'point';     % apex pencil (Up=Dn=tip)

    % ---- Misc ----
    C.eps1     = 1e-8;
    C.fact     = 30;             % visualization scale
    C.plotGeom = false;
    C.plotMesh = false;
    C.stressDPI = 600;

    % ----------------- apply user overrides -----------------
    if ~isempty(varargin)
        C = apply_overrides(C, varargin{:});
        % ensure nLeg consistent if user overwrote L/theta
        C.nLeg = numel(C.L);
    end

    % ----------------- build polyline vertices -----------------
    C.Pmid = build_polyline(C.P0, C.L, C.theta);
    
    % --------% initial trial cohesive length (baseline probe leg) -------
    C.delta0 = norm(C.Pmid(end,:) - C.Pmid(end-1,:)); 

    % Local frame at the LAST segment (used by your current Dt-mouth extraction)
    C.Q = local_frame_last_segment(C.Pmid);

    % ---- Pencil half-width (channel half-width) ----
    % Use LAST leg length as cohesive leg length proxy:
    Llast = C.L(end);
    hmin_last = Llast / C.ncoh;
    C.chw = hmin_last / 16;   % small relative to local feature size

    % ----------------- material -----------------
    C.nu = 0.3;
    C.E2 = 4e3;
    C.Dmat = C.E2*[[1-C.nu, C.nu, 0];
                   [C.nu, 1-C.nu, 0];
                   [0, 0, (1-2*C.nu)/2]]/(1+C.nu)/(1-2*C.nu);
    C.G12 = C.E2/2/(1+C.nu);

    % ----------------- TSL params -----------------
    C.phi    = [800, 1200];
    C.rphi   = 0.4;
    C.a1     = 0.002;
    C.a2     = [.7, .5];
    C.rEs    = 200;
    C.sigmax = [C.E2, C.G12]/C.rEs;

    % area under normalized traction curve (vectorized)
    c_area   = (3 - 2*C.a1 + 3*C.a2)/6;
    C.Dmax   = C.phi ./ c_area ./ (C.sigmax*1e6);

    % load guess
    C.sig_guess = C.sigmax(1)/8;

    % Convenience: store the swept angle explicitly for logging
    C.theta_k = C.theta(end);
end

% ========================================================================

function P = build_polyline(P0, L, theta)
%BUILD_POLYLINE vertices from mouth point + (length, angle) legs
    nLeg = numel(L);
    P = zeros(nLeg+1, 2);
    P(1,:) = P0;
    for j = 1:nLeg
        P(j+1,:) = P(j,:) + L(j)*[cos(theta(j)), sin(theta(j))];
    end
end

function Q = local_frame_last_segment(Pmid)
%LOCAL_FRAME_LAST_SEGMENT local basis from the last polyline segment
% Returns Q such that [d_local] = Q * [jump_global]
% (keep this consistent with your existing post-processing convention).
    p1 = Pmid(end-1,:);
    p2 = Pmid(end,:);
    v  = p2 - p1;

    nv = norm(v);
    if nv < eps
        error('local_frame_last_segment: last segment has near-zero length.');
    end

    t = v / nv;               % tangent
    n = [-t(2), t(1)];         % normal (90 deg CCW)

    Q = [t; n];
end

function C = apply_overrides(C, varargin)
%APPLY_OVERRIDES name-value parser for common overrides
    if mod(numel(varargin),2) ~= 0
        error('cfg_min_geom: overrides must be name-value pairs.');
    end
    for k = 1:2:numel(varargin)
        name = varargin{k};
        val  = varargin{k+1};
        if ~ischar(name) && ~isstring(name)
            error('cfg_min_geom: override names must be char/string.');
        end
        name = char(lower(string(name)));

        switch name
            case 'p0'
                C.P0 = val;
            case 'l'
                C.L = val(:);
            case 'theta'
                C.theta = val(:);
            case 'ncoh'
                C.ncoh = val;
            case 'hmax_ratio'
                C.hmax_ratio = val;
            case 'hgrad'
                C.hgrad = val;
            case 'a'
                C.A = val;
            case 'b'
                C.B = val;
            case 'join'
                C.join = val;
            case 'miter_limit'
                C.miter_limit = val;
            case 'corner_tol'
                C.corner_tol = val;
            case 'tip'
                C.tip = val;
            case 'plotgeom'
                C.plotGeom = logical(val);
            case 'plotmesh'
                C.plotMesh = logical(val);
            case 'fact'
                C.fact = val;
            otherwise
                % allow arbitrary fields without breaking (handy during development)
                C.(varargin{k}) = val;
        end
    end

    % sanity if user overwrote L/theta
    if numel(C.L) ~= numel(C.theta)
        error('cfg_min_geom: numel(L) must equal numel(theta).');
    end
end
