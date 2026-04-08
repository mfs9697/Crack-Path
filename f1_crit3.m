function [F_aug, J_aug] = f1_crit3(x, mesh, cz, Stif, elod, cohes)
%F1_CRIT3  Critical-state augmented equilibrium for conservative CZM
%          on a polyline crack with a cohesive segment.
%
% Unknowns: x = [u; sig], u in R^{ndof}, sig scalar.
% Constraint: Dn at mouth equals Dmax(1).
%
% Cohesive traction is applied on ALL cohesive quadratic elements
% (including the last element that contains the crack tip), while enforcing
% a traction-free crack tip by setting the tip-station traction AND its
% consistent tangent to zero.
%
% Inputs:
%   mesh.coord   : [N x 2] node coordinates
%   mesh.Pmid0   : [Np x 2] crack midline polyline; last segment is cohesive leg
%   cz           : cohesive-law struct for tsl2
%                 Required fields:
%                   cz.Dmax   [Dnmax Dtmax]
%                   cz.sigmax [sig_nmax sig_tmax]
%                   + law parameters used by tsl2
%   Stif         : bulk stiffness (ndof x ndof sparse)
%   elod         : [nload x 2] = [node_id, weight] for vertical DOF loads
%   cohes        : (2*ncoh+1) x 2 station ids [upper lower], mouth->tip
%                 (interleaved vertex, mid, vertex, ...)
%
% Outputs:
%   F_aug  : [ndof+1 x 1] residual = [equilibrium; constraint]
%   J_aug  : [(ndof+1) x (ndof+1)] Jacobian

if nargin < 6 || isempty(cohes)
    if isfield(cz,'cohes')
        cohes = cz.cohes;
    else
        error('f1_crit3:MissingCohes', 'Provide cohes (or set cz.cohes).');
    end
end

% ---- construct local frame m = [t; n] from the cohesive leg (last polyline segment)
m=cz.Q;

coord = mesh.coord;
ndof  = 2*size(coord,1);
u     = x(1:ndof);
sig   = x(end);

Dmax  = cz.Dmax(:).';      % [Dnmax Dtmax]
sigmx = cz.sigmax(:).';    % [sig_nmax sig_t_max]

% ---- external unit load vector (uy DOFs only)
F_ext = sparse(2*elod(:,1), 1, elod(:,2), ndof, 1);

% ---- quadratic edge “mass” matrix (3-node quadratic line)
R  = [4 2 -1; 2 16 2; -1 2 4]/30;
I3 = eye(3);

% ---- cohesive discretization sizes
mcoh      = size(cohes,1);      % number of stations (interleaved v,m,v,...)
ncoh_elem = (mcoh - 1)/2;       % number of quadratic cohesive elements

% ---- per-station storage
sigma_local     = zeros(mcoh, 2);    % [t, n] physical traction at stations
dsigma_dD_local = cell(mcoh, 1);     % [t;n] x [Dt;Dn] physical tangents
for iS = 1:mcoh
    dsigma_dD_local{iS} = zeros(2,2);
end

% =========================
% 1) Station-wise cohesive law evaluation (ALL stations)
% =========================
for iS = 1:mcoh
    iu = cohes(iS,1);
    il = cohes(iS,2);

    % global jump (top - bottom)
    Dg = u([2*iu-1, 2*iu]) - u([2*il-1, 2*il]);

    % local jump [Dt; Dn]
    Dl = m * Dg;

    xi_n_signed = Dl(2) / Dmax(1);
    xi_t_signed = Dl(1) / Dmax(2);

    % strict sign for Dt (nondifferentiable at 0; keep as in current project)
    sgn_t = sign(xi_t_signed);

    % conservative law uses abs(xi_t) internally; returns magnitudes:
    %   f_nt = [fn; ft_mag], M_nt is symmetric in [n,t] ordering
    [f_nt, M_nt] = tsl2([xi_n_signed, xi_t_signed], cz);

    % physical tractions in [t, n] order for assembly
    ft = (f_nt(2) * sgn_t) * sigmx(2);
    fn =  f_nt(1)           * sigmx(1);
    sigma_local(iS,:) = [ft, fn];

    % physical tangent for magnitudes: d[sig_n; |sig_t|]/d[Dn; |Dt|] in [n,t]
    Sphys = diag([ sqrt(sigmx(1)/Dmax(1)),  sqrt(sigmx(2)/Dmax(2)) ]); % [n,t]
    K_mag = Sphys * M_nt * Sphys;  % [n,t] x [Dn,|Dt|]

    % map |Dt| -> Dt_signed: multiply Dt-column by sgn_t
    K_mag(:,2) = K_mag(:,2) * sgn_t;

    % map sigma_t = sgn_t * |sigma_t|: multiply shear-equation row by sgn_t
    K_mag(2,:) = K_mag(2,:) * sgn_t;

    % reorder to [t;n] x [Dt;Dn] for assembly
    dsigma_dD_local{iS} = K_mag([2,1],[2,1]);
end

% =========================
% 1b) Enforce traction-free crack tip (station-level, consistent)
% =========================
itip = mcoh;  % last station in mouth->tip ordering is the crack tip station
sigma_local(itip,:)      = [0, 0];
dsigma_dD_local{itip}    = zeros(2,2);

% =========================
% 2) Assemble cohesive residual and Jacobian (ALL elements)
% =========================
J_coh = sparse(ndof, ndof);
clod  = zeros(12*ncoh_elem, 2);   % [dof_id, value]
im3   = 1;

for iC = 1:ncoh_elem
    ind = 2*(iC-1) + (1:3);       % station rows [v, m, v] for element iC

    up  = cohes(ind,1);
    lo  = flipud(cohes(ind,2));   % reverse lower to ensure consistent CCW ordering

    % element length (top end nodes)
    le       = norm(coord(up(3),:) - coord(up(1),:));
    R_scaled = le * R;

    % local->global mapping for nodal forces
    G = kron(R_scaled, m.');      % 6x6

    % residual on upper face (minus sign because internal forces)
    s_up = -reshape(sigma_local(ind,:).', [], 1);   % [t1;n1;t2;n2;t3;n3]
    t_up = G * s_up;                                 % 6x1
    t_lo = -t_up([5,6,3,4,1,2]);                     % action-reaction, reordered

    top_dofs = reshape([2*up-1, 2*up].', [], 1);
    bot_dofs = reshape([2*lo-1, 2*lo].', [], 1);

    clod(im3:im3+11,:) = [top_dofs, t_up; bot_dofs, t_lo];
    im3 = im3 + 12;

    % element Jacobian
    Dblk = blkdiag(dsigma_dD_local{ind(1)}, dsigma_dD_local{ind(2)}, dsigma_dD_local{ind(3)});
    m_e  = kron(I3, m);           % 6x6 mapping global jump -> local blocks

    % jump selector b: maps [u_top; u_bot] -> stacked local jumps at 3 stations
    B1 = [ eye(2),             zeros(2,8),    -eye(2)          ];
    B2 = [ zeros(2),           eye(2), zeros(2,4), -eye(2),    zeros(2) ];
    B3 = [ zeros(2,4),         eye(2),         -eye(2),        zeros(2,4) ];
    b  = [B1; B2; B3];           % 6x12

    J_up = -G * Dblk * m_e * b;
    J_lo = -J_up([5,6,3,4,1,2],:);

    dofs = [top_dofs; bot_dofs];
    J_e  = [J_up; J_lo];

    J_coh(dofs,dofs) = J_coh(dofs,dofs) + J_e;
end

F_coh = sparse(clod(:,1), 1, clod(:,2), ndof, 1);

% =========================
% 3) Bulk equilibrium + homotopy scaling
% =========================
F_eq = Stif*u - (sig*F_ext + F_coh);
J_u  = Stif   - J_coh;

% =========================
% 4) Mouth constraint: Dn_mouth = Dmax(1)
% =========================
iu   = cohes(1,1); il = cohes(1,2);   % mouth station
nvec = m(2,:);                        % local normal

g = nvec * ( u([2*iu-1;2*iu]) - u([2*il-1;2*il]) ) - Dmax(1);

F_aug = [F_eq; g];
J_aug = [J_u,                               -F_ext; ...
         sparse(1, [2*iu-1,2*iu,2*il-1,2*il], [nvec, -nvec], 1, ndof),  0];
end
