function [coord1, str] = StressExt(mesh, U, mat, quad, fact)
%STRESSEXT  Recover nodal stresses for T6 mesh by GP evaluation and
%           extrapolation + area-weighted averaging.
%
% Robust fixes vs older version:
%   - use NaN occupancy instead of "==0"
%   - abs(triangle area) for weights
%   - overflow guard on per-node adjacency capacity

coord0    = mesh.coord;
connect0  = mesh.connect;
Dmat      = mat.D;

nip2  = quad.nip2;
xip2  = quad.xip2;
Nextr = quad.Nextr;

nnod  = size(coord0, 1);
nelem = size(connect0, 1);

coord1 = coord0 + reshape(U, 2, [])' * fact;

% ---- guards
if size(Nextr,1) ~= 6 || size(Nextr,2) ~= nip2
    error('StressExt:BadNextr', 'Expected Nextr to be 6 x nip2.');
end
if ~isequal(size(xip2), [2, nip2])
    error('StressExt:BadXip2', 'Expected xip2 to be 2 x nip2.');
end
if size(connect0,2) ~= 6
    error('StressExt:NotT6', 'mesh.connect must be nelem x 6 for T6.');
end

% ---- storage
str      = zeros(nnod, 3);        % final nodal stresses
str_ip   = zeros(nip2, nelem, 3); % GP stresses
aelem    = zeros(nelem, 1);       % triangle areas
eldof    = zeros(12, 1);          % element dofs (2*6)

% adjacency capacity (increase if needed)
CAP = 30;
str_adj  = NaN(nnod, 3, CAP);     % per-node stress contributions (NaN = empty)
elmn_adj = zeros(nnod, CAP);      % element ids per contribution

for elem = 1:nelem
    nodes = connect0(elem,:);

    % element dofs
    for n = 1:6
        eldof(2*n-1:2*n) = [2*nodes(n)-1; 2*nodes(n)];
    end
    u_elem = U(eldof);

    % Gauss-point stresses
    for ig = 1:nip2
        [B2, ~] = BN_local(xip2(:,ig), coord0(nodes,:));
        strain = B2 * u_elem;          % [exx; eyy; gxy]
        stress = Dmat * strain;        % [sxx; syy; sxy]
        str_ip(ig, elem, :) = stress;
    end

    % Extrapolate each component to the 6 T6 nodes and store per-node contributions
    for s = 1:3
        stress_gp    = str_ip(:, elem, s);
        stress_nodal = Nextr * stress_gp;  % 6x1

        for a = 1:6
            node = nodes(a);

            slot = find(isnan(str_adj(node, s, :)), 1, 'first');
            if isempty(slot)
                error('StressExt:AdjOverflow', ...
                    'Adjacency overflow at node %d. Increase CAP.', node);
            end

            str_adj(node, s, slot) = stress_nodal(a);
            elmn_adj(node, slot)   = elem;
        end
    end

    % Triangle area from 3 vertex nodes (use abs for safety)
    x = coord0(nodes(1:3), :)';
    aelem(elem) = abs(det([x; 1 1 1])) / 2;
end

% Area-weighted nodal averaging
for n = 1:nnod
    for s = 1:3
        vals = squeeze(str_adj(n, s, :));
        k    = find(~isnan(vals));
        if isempty(k), continue; end

        elems = elmn_adj(n, k);
        areas = aelem(elems);

        % guard against zero total area (should not happen)
        denom = sum(areas);
        if denom <= 0
            continue;
        end

        str(n, s) = (vals(k)' * areas) / denom;
    end
end
end

% ----- local helper -----
function [B,Det] = BN_local(xi0, X)
% T6 B-matrix
xi  = [xi0; 1 - sum(xi0)];
Nap = [ 4*xi(1)-1, 0,          1-4*xi(3), 4*xi(2),        -4*xi(2),       4*xi(3)-4*xi(1);
        0,          4*xi(2)-1, 1-4*xi(3), 4*xi(1),         4*xi(3)-4*xi(2), -4*xi(1) ];

dxdxi = Nap * X; Det = det(dxdxi);
N1    = [ dxdxi(2,2), -dxdxi(1,2);
         -dxdxi(2,1),  dxdxi(1,1) ] / Det * Nap;

eldf  = 12; inx=(2:2:eldf)'-1; iny=inx+1;
B = zeros(3, eldf);
B(1, inx) = N1(1,:);
B(2, iny) = N1(2,:);
B(3, inx) = N1(2,:);
B(3, iny) = N1(1,:);
end
