function [mesh, mat, quad, elod, Stif, fixvar, id_tip] = geom_pencil_1leg(C,arrow)
%GEOM_PENCIL_1LEG  Geometry + mesh for a 1-leg straight crack (LEFM)
% with a pencil-like thin channel included in the PDE boundary.
%
% Channel:
%   - from p0 (mouth) to p1 (1-arrow a) : constant half-width = C.chw
%   - from p1 to p2 (tip)            : faces collapse to a single tip point
%
% Inputs in struct C:
%   A, B      - strip half-widths: domain [0, A] x [-B, B]
%   a         - crack length
%   alf1      - crack angle (radians) from +x axis
%   chw       - half-width of thin channel around crack
%   hmax      - target max element size away from crack
%   hgrad     - mesh gradation parameter
%   eps1      - smoothing parameter for edge loads
%   E2, nu    - elastic constants (plane strain)
%   G12       - shear modulus (if needed elsewhere)
%   Dmat      - 3x3 elasticity matrix for plane strain (global Dmat)
%
% Outputs:
%   V      - [3x2] crack points [p0; p1; p2] (mouth, 0.99a, tip)
%   mesh   - struct with T6 and parent T3 meshes
%   mat    - material struct
%   quad   - quadrature struct (nip2, xip2, w2, Nextr)
%   elod   - nodal edge loads (top/bottom)
%   Stif   - global stiffness matrix (BCs applied)
%   fixvar - fixed DOF indices

here = fileparts(mfilename('fullpath'));
utilspath = fullfile(here, 'utils');
if exist(fullfile(utilspath,'distanceFromLine.m'),'file') && isempty(which('distanceFromLine'))
    addpath(utilspath);
end

% ---------- inputs ----------
A    = C.A;
B    = C.B;
a    = C.a;
chw  = C.chw;
hmax = C.hmax;
hgrad= C.hgrad;
eps1 = C.eps1;

% ---------- crack points ----------
p0 = [0, 0];                          % mouth
p2 = a * [cos(C.alf1), sin(C.alf1)];  % tip
AB = p2 - p0;
L  = max(1e-30, norm(AB));
e  = AB / L;                          % unit along crack

s_knee = (1-arrow) * L;
p1     = p0 + s_knee * e;             % "knee" at 0.99 a

% ---------- Step 1: build polygon with pencil channel + outer rectangle ----------
[Bound, G] = build_polygon_pencil_1leg(p0, p1, p2, A, B, chw);

% CSG: a single polygon region (same style as original geom_pencil.m)
Gcol = [2; size(Bound,1); Bound(:,1); Bound(:,2)];
g    = decsg(Gcol);
mdl  = createpde();
geometryFromEdges(mdl, g);

% ---------- Step 2: T3 mesh with tip refinement (hmin = a*arrow) ----------
hmin  = max(a*arrow, 1e-3);            
hmaxL = max(hmax, 2*hmin);
hgradL= max(1.05, min(hgrad, 1.7));

% refine near the geometry vertex closest to the crack tip
verts = mdl.Geometry.Vertices;        % 2 x Nvert
if size(verts,1) == 2 && size(verts,2) ~= 2
    verts = verts.';                  % Nvert x 2
elseif size(verts,2) == 2 && size(verts,1) ~= 2
    % already Nvert x 2
else
    error('Unexpected size of mdl.Geometry.Vertices.');
end
[~, icrt] = min(sum((verts - p2).^2, 2));

Mesh = generateMesh(mdl, ...
    'Hmin',     hmin, ...
    'Hmax',     hmaxL, ...
    'Hgrad',    hgradL, ...
    'Hvertex', {icrt, hmin/10}, ...
    'GeometricOrder', 'linear');
%figure(3); clf(3); pdemesh(mdl)

coord3   = Mesh.Nodes.';      % T3 coords (pre-collapse)
connect3 = Mesh.Elements.';   % T3 connectivity

% ---------- Step 3: find channel faces on T3 via segments ----------
%
% Channel boundary in terms of G:
%   upper: up_mouth -> up_knee -> tip
%   lower: tip      -> dn_knee -> dn_mouth
%

leg_up = [];
leg_dn = [];

% upper: p0 -> p1 -> p2
Upts = [G.up_mouth; G.up_knee; G.p2];
for k = 1:(size(Upts,1)-1)
    iu = face_ids_on_segment(coord3, Upts(k,:), Upts(k+1,:));
    if k > 1 && ~isempty(iu)
        iu = iu(2:end);  % drop duplicate at the joint
    end
    leg_up = [leg_up; iu]; %#ok<AGROW>
end

% lower: p2 -> p1 -> p0
Dpts = [G.p2; G.dn_knee; G.dn_mouth];
for k = 1:(size(Dpts,1)-1)
    id = face_ids_on_segment(coord3, Dpts(k,:), Dpts(k+1,:));
    if k > 1 && ~isempty(id)
        id = id(2:end);  % drop duplicate
    end
    leg_dn = [leg_dn; id]; %#ok<AGROW>
end

% ---------- Step 4: collapse channel nodes to pencil-shaped crack line ----------
coord_lin = coord3;

rel_all = coord_lin - p0;
s_all   = rel_all * e.';       % along-crack coordinate (0..L)

ids_chan = unique([leg_up; leg_dn]);  % all channel-face nodes
id_tip = intersect(leg_up, leg_dn);
if ~isempty(ids_chan)
    s_nodes = s_all(ids_chan);
    proj    = coord_lin(ids_chan,:);  % initialize with current positions

    % --- classify tip vs midline by closeness to L ---
    Ltol   = 1e-6 * max(1.0, L);      % small absolute tolerance
    isTip  = (s_nodes >= L - Ltol);   % ONLY the true tip node(s)
    isMid  = ~isTip;

    % midline projection for all non-tip channel nodes
    s_mid = max(0, min(L, s_nodes(isMid)));
    proj(isMid,:) = p0 + s_mid .* e;

    % tip: project exactly to p2
    proj(isTip,:) = repmat(p2, sum(isTip), 1);

    % write back
    coord_lin(ids_chan,:) = proj;
end

%figure(4); clf(4); triplot(connect3, coord_lin(:,1), coord_lin(:,2)); hold on
%axis equal

% ---------- Step 5: clean triangles + orient CCW ----------
Atri = triAreasSigned(connect3, coord_lin);
keep = abs(Atri) > 1e-15;
connect3 = connect3(keep,:); Atri = Atri(keep);
cw = Atri < 0;
if any(cw)
    tmp = connect3(cw,2);
    connect3(cw,2) = connect3(cw,3);
    connect3(cw,3) = tmp;
end

% ---------- Step 6: upgrade T3 -> T6 ----------
coord   = coord_lin;
connect = connect3;
[coord, connect] = T3toT6(coord, connect);

% ---------- pack mesh struct ----------
mesh = struct();
mesh.coord    = coord;      % T6
mesh.connect  = connect;    % T6
mesh.coord3   = coord_lin;  % T3 (collapsed)
mesh.connect3 = connect3;   % T3

% ---------- material ----------
mat = struct('E',   C.E2, ...
             'nu',  C.nu, ...
             'G12', C.G12, ...
             'D',   C.Dmat);

% ---------- quadrature ----------
[nip2, xip2, w2, Nextr] = integr();
quad = struct('nip2', nip2, 'xip2', xip2, 'w2', w2, 'Nextr', Nextr);

% ---------- edge loads on top/bottom for unit nominal stress ----------
elod = edge_loads_T6(coord, B, eps1);

% ---------- essential BCs (same as in your 2-leg code) ----------
fix_pts = [A, 0; 0, B; 0, -B];
fix = zeros(3,1);
for i = 1:3
    [~,fix(i)] = min((coord(:,1)-fix_pts(i,1)).^2 + (coord(:,2)-fix_pts(i,2)).^2);
end
fixvar = [2*fix(1); 2*fix(2)-1; 2*fix(3)-1];

% ---------- assemble stiffness ----------
Stif = stif_assem(mesh, mat, quad, fixvar);

% --- diagnostics: look for completely free DOFs ---
rowSum = sum(abs(Stif), 2);
zeroRows = find(rowSum == 0);

if ~isempty(zeroRows)
    fprintf('WARNING: %d zero rows in Stif (orphan DOFs)\n', numel(zeroRows));
end


% ---------- optional plots ----------
if ~isfield(C,'plotMesh') || C.plotMesh
    figure(3); clf;

    % (a) pre-collapse: polygon + T3 mesh
    subplot(1,2,1)
    triplot(connect3, coord3(:,1), coord3(:,2)); hold on
    Cpoly = [Bound; Bound(1,:)];
    plot(Cpoly(:,1), Cpoly(:,2), 'r-', 'LineWidth', 1.0)
    box on; title('(a) Pencil channel in polygon (T3)')
    axis equal; xlim([-0.1*A, 1.1*A]); ylim([-1.1*B, 1.1*B]);

    % (b) post-collapse: crack line visible, tip is a single point
    subplot(1,2,2)
    triplot(connect(:,1:3), coord(:,1), coord(:,2)); hold on
    rect = [0 -B; A -B; A B; 0 B; 0 -B];
    plot(rect(:,1), rect(:,2), 'r-', 'LineWidth', 1.0)
    plot([p0(1) p2(1)], [p0(2) p2(2)], 'r-', 'LineWidth', 2)
    box on; title('(b) Crack mesh (T6, pencil tip)')
    axis equal; xlim([-0.1*A, 1.1*A]); ylim([-1.1*B, 1.1*B]);
end

end

% ========================================================================
%                               HELPERS
% ========================================================================

function [Bound, G] = build_polygon_pencil_1leg(p0, p1, p2, A, B, w)
%BUILD_POLYGON_PENCIL_1LEG  Polygon with pencil channel + rectangle.
%
% upper chain: up_mouth -> up_knee -> tip
% lower chain: tip -> dn_knee -> dn_mouth
% then outer rectangle [0,A] x [-B,B].

e = p2 - p0; 
e = e / max(norm(e), eps);
n = [-e(2), e(1)];   % left normal

up_mouth = p0 + w*n;
dn_mouth = p0 - w*n;
up_knee  = p1 + w*n;
dn_knee  = p1 - w*n;
tip      = p2;

upper_chain = [up_mouth; up_knee; tip];
lower_chain = [dn_knee; dn_mouth];

Bound = [ upper_chain;
          lower_chain;
          0, -B;
          A, -B;
          A,  B;
          0,  B ];

if signed_area(Bound) < 0
    Bound = flipud(Bound);
end

G.up_mouth = up_mouth;
G.dn_mouth = dn_mouth;
G.up_knee  = up_knee;
G.dn_knee  = dn_knee;
G.p2       = tip;
end

function ids = face_ids_on_segment(P, A, B)
%FACE_IDS_ON_SEGMENT: node indices in P that lie on segment A->B.

AB = B - A;
L2 = max(1e-32, dot(AB,AB));

if exist('distanceFromLine','file') == 2
    d = distanceFromLine(A,B,P);
else
    n = [-AB(2), AB(1)];
    nn = norm(n);
    if nn == 0
        d = zeros(size(P,1),1);
    else
        d = (P - A) * (n.'/nn);
    end
end

t = ((P - A) * AB.') / L2;
mask = (abs(d) < 1e-12) & (t >= -1e-10) & (t <= 1+1e-10);
ids  = find(mask);
if isempty(ids), return, end

s = vecnorm(P(ids,:) - A, 2, 2);
ids = sortrows([ids, s], 2);
ids = ids(:,1);
end

function A = signed_area(P)
x = P(:,1); y = P(:,2);
A = 0.5*sum(x.*y([2:end 1]) - x([2:end 1]).*y);
end

function A = triAreasSigned(T, X)
v1 = X(T(:,1),:); v2 = X(T(:,2),:); v3 = X(T(:,3),:);
A  = 0.5*((v2(:,1)-v1(:,1)).*(v3(:,2)-v1(:,2)) - ...
          (v2(:,2)-v1(:,2)).*(v3(:,1)-v1(:,1)));
end
