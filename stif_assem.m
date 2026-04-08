function Stif = stif_assem(mesh, mat, quad, fixvar)
% STIF_ASSEM  Assemble global stiffness for T6 triangles (nelvert=6).
% Inputs
%   mesh   : struct with fields
%              - coord  [N x 2]
%              - connect[Ne x 6]
%   mat    : struct with fields
%              - D  [3 x 3] plane strain/stress matrix
%   quad   : struct with fields
%              - nip2, xip2 [2 x nip2], w2 [1 x nip2]
%   fixvar : column vector of global DOF indices to clamp (rows -> I)
%
% Output
%   Stif   : sparse global stiffness (Dirichlet rows set to identity)

% ---- unpack ----
coord   = mesh.coord;
connect = mesh.connect;
Dmat    = mat.D;

nip2 = quad.nip2;
xip2 = quad.xip2;
w2   = quad.w2;

% ---- dimensions & prealloc ----
nelvert = 6; eldf = 2*nelvert;
ndof = 2*size(coord,1);

% heuristic for triplet preallocation (10% of dense; adjust if needed)
m0 = max(1, round(0.05*ndof^2));
rw = zeros(m0,1); cl = rw; st = rw; k = 0;

% ---- element loop ----
for e = 1:size(connect,1)
    nodes = connect(e,:);
    X     = coord(nodes,:);    % [6 x 2]

    % element stiffness
    Ke = zeros(eldf);
    for j = 1:nip2
        [B, DetJ] = BN_local(xip2(:,j), X);   % local helper below
        % (Optional) guard against inverted elements:
        if DetJ <= 0, error('Inverted or degenerate element at e=%d', e); end
        Ke = Ke + w2(j) * (B.' * Dmat * B) * (DetJ/2);
    end

    % assemble into triplets
    % map local (i,alpha) to global dof: dof = 2*node-1+(alpha-1), alpha in {1,2}
    edofs = zeros(eldf,1);
    for n = 1:nelvert
        edofs(2*n-1:2*n) = [2*nodes(n)-1; 2*nodes(n)];
    end

    % add all entries of Ke
    [ii, jj] = ndgrid(1:eldf, 1:eldf);
    nn = numel(ii);
    if k + nn > numel(rw)
        % grow storage (x1.8)
        grow = max(nn, round(0.8*numel(rw)));
        rw = [rw; zeros(grow,1)];
        cl = [cl; zeros(grow,1)];
        st = [st; zeros(grow,1)];
    end
    rw(k+(1:nn)) = edofs(ii(:));
    cl(k+(1:nn)) = edofs(jj(:));
    st(k+(1:nn)) = Ke(sub2ind([eldf,eldf], ii(:), jj(:)));
    k = k + nn;
end

% ---- finalize sparse ----
rw = rw(1:k); cl = cl(1:k); st = st(1:k);
Stif = sparse(rw, cl, st, ndof, ndof);

% ---- impose Dirichlet rows as identity (standard row clamping) ----
if ~isempty(fixvar)
    % zero out existing entries in those rows
    Stif(fixvar, :) = 0;
    % set diagonal to 1 on clamped dofs
    Stif = Stif + sparse(fixvar, fixvar, ones(numel(fixvar),1), ndof, ndof);
end

end

