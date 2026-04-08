function DispStressExt(mesh, U, cz, mat, quad, nstress, fact)
% DispStressExt
% Visualizes:
%  - stress field (normalized component or von Mises)
%  - normalized openings and cohesive tractions on both faces
%  - TSL curve overlay (using the same sign convention as solver)
%
% Convention (locked):
%   cz.cohes(1,:)   = CZ mouth station (start of cohesive leg)
%   cz.cohes(end,:) = tip station
%   xlocal measured from CZ mouth along the cohesive-leg tangent
%   traction-free tip enforced in the plotted TSL tractions

% --- required inputs ---
Dmat = mat.D; E = mat.E; nu = mat.nu;

% --- Stress field ---
[coord1, sig_node] = StressExt(mesh, U, mat, quad, fact);
coord   = mesh.coord;
connect = mesh.connect;

nelem = size(connect,1);
vert  = [1,4,2,5,3,6];   % node order for plotting T6

% normalization for color field
if     nstress == 3
    smax = max(1e-12, cz.sigmax(2));      % note: used for shear component plot
elseif nstress == 4
    smax = 1;                             % von Mises absolute
else
    smax = max(1e-12, cz.sigmax(1));      % normal components
end

X = zeros(nelem, numel(vert));
Y = X; Z = X;

for e = 1:nelem
    sig0 = sig_node(connect(e,vert),:);  % [sxx, syy, sxy]

    if nstress < 4
        Z(e,:) = sig0(:,nstress) / smax;
    else
        % 3D-equivalent von Mises (plane strain sigma_z)
        eps0 = (Dmat \ sig0')';    % strains
        sigx = sig0(:,1); sigy = sig0(:,2); sigxy = sig0(:,3);
        sigz = E*nu/(1+nu)/(1-2*nu) * (eps0(:,1) + eps0(:,2));
        vm   = sqrt((sigx - sigy).^2 + (sigy - sigz).^2 + (sigz - sigx).^2 + 6*sigxy.^2)/sqrt(2);
        Z(e,:) = vm;
    end

    X(e,:) = coord1(connect(e,vert),1);
    Y(e,:) = coord1(connect(e,vert),2);
end

figure(2); clf(2)
ax = gca; hold(ax,'on'); view(ax,2);
set(ax,'DataAspectRatio',[1 1 1]);

col = linspace(0,1,21)'; cmap = jet(numel(col)-1);
clim([col(1) col(end)]); colormap(cmap); colorbar('Ticks',col);
patch(X',Y',Z','EdgeColor',.6*[1,1,1],'LineWidth',.05,'FaceColor','interp');

if     nstress == 1
    title('$\sigma_{xx} / \sigma_{(n)}^{\max}$','Interpreter','latex');
elseif nstress == 2
    title('$\sigma_{yy} / \sigma_{(n)}^{\max}$','Interpreter','latex');
elseif nstress == 3
    title('$\sigma_{xy} / \sigma_{(t)}^{\max}$','Interpreter','latex');
else
    title('von Mises (absolute)');
end
box on

% ------------------ CZ stations & windowing ------------------
cohes = cz.cohes;                 % (2*ncoh+1) x 2, mouth->tip
mcoh  = size(cohes,1);

% Local frame Q = [t; n] from last polyline segment (consistent with solver)
if ~isfield(mesh,'Pmid0') || size(mesh.Pmid0,1) < 2
    error('DispStressExt:MissingPmid0', 'mesh.Pmid0 is required to build the CZ frame.');
end
p1 = mesh.Pmid0(end-1,:);
p2 = mesh.Pmid0(end,:);
tv = p2 - p1;
nv = norm(tv);
if nv < 1e-14
    error('DispStressExt:DegenerateLeg', 'Last segment of mesh.Pmid0 has zero length.');
end
t = tv / nv;
n = [-t(2), t(1)];
Q = [t; n];                       % rows: [t; n]

% CZ mouth and tip station ids
iu0 = cohes(1,1);    il0 = cohes(1,2);      % CZ mouth
iu1 = cohes(end,1);  il1 = cohes(end,2);    % tip

p_mouth0 = 0.5*(coord(iu0,:) + coord(il0,:)); % CZ mouth midpoint (undeformed)
p_tip0   = 0.5*(coord(iu1,:) + coord(il1,:)); % tip midpoint (undeformed)
da_est   = norm(p_tip0 - p_mouth0);           % cohesive-leg length estimate

% center view on the DEFORMED tip
pu = coord(iu1,:) + fact * [U(2*iu1-1), U(2*iu1)];
pl = coord(il1,:) + fact * [U(2*il1-1), U(2*il1)];
p_tip = 0.5*(pu + pl);

indent = [2 1.1] * da_est;
xlim(ax, [p_tip(1)-indent(1), p_tip(1)+indent(2)]);
ylim(ax, [p_tip(2)-.6*indent(1), p_tip(2)+.6*indent(1)]);
box(ax,'on');

xlabel('$x$, m','Interpreter','latex');
ylabel('$y$, m','Interpreter','latex');

% ------------------ Cohesive openings & tractions (ALL stations) ------------------
xlocal = zeros(mcoh,1);
Dn     = zeros(mcoh,1);   % normalized opening (signed in compression/opening convention)
Dt     = zeros(mcoh,1);   % normalized tangential jump (signed)

% From bulk stresses (normalized): "bulk-projected" traction (qualitative)
Tn_top_bulk = zeros(mcoh,1);  Tt_top_bulk = zeros(mcoh,1);
Tn_bot_bulk = zeros(mcoh,1);  Tt_bot_bulk = zeros(mcoh,1);

% From TSL via jump (normalized)
Tn_top_tsl  = zeros(mcoh,1);  Tt_top_tsl  = zeros(mcoh,1);
Tn_bot_tsl  = zeros(mcoh,1);  Tt_bot_tsl  = zeros(mcoh,1);

Dnmax = max(1e-12, cz.Dmax(1));
Dtmax = max(1e-12, cz.Dmax(2));
snmax = max(1e-12, cz.sigmax(1));
stmax = max(1e-12, cz.sigmax(2));

for iS = 1:mcoh
    iu = cohes(iS,1);
    il = cohes(iS,2);

    % local abscissa from CZ mouth (undeformed midpoint)
    pmid      = 0.5*(coord(iu,:) + coord(il,:));
    xlocal(iS)= Q(1,:) * (pmid - p_mouth0).';

    % openings (normalized) from displacement jump (local)
    Uu   = U([2*iu-1, 2*iu]);
    Ul   = U([2*il-1, 2*il]);
    Dloc = Q * (Uu - Ul);               % [Dt; Dn] in local frame

    Dt(iS) = Dloc(1) / Dtmax;
    Dn(iS) = Dloc(2) / Dnmax;

    % bulk-projected traction using recovered nodal stress (qualitative)
    st  = sig_node(iu,:).';
    Sg  = [st(1) st(3); st(3) st(2)];
    Sl  = Q * Sg * Q.';
    t_top_loc = -Sl * [0;1];            % traction on +n
    Tt_top_bulk(iS) = t_top_loc(1) / stmax;
    Tn_top_bulk(iS) = t_top_loc(2) / snmax;

    sb  = sig_node(il,:).';
    Sgb = [sb(1) sb(3); sb(3) sb(2)];
    Slb = Q * Sgb * Q.';
    t_bot_loc =  Slb * [0;1];
    Tt_bot_bulk(iS) = t_bot_loc(1) / stmax;
    Tn_bot_bulk(iS) = t_bot_loc(2) / snmax;

    % TSL tractions from jump (normalized) using solver-like shear sign handling
    xi_n  = Dn(iS);
    xi_t  = Dt(iS);
    sgn_t = sign(xi_t);

    % Expect tsl2 to return magnitudes: f_nt = [fn; ft_mag] (dimensionless)
    [f_nt, ~] = tsl2([xi_n, xi_t], cz);

    fn = f_nt(1);
    ft = f_nt(2) * sgn_t;

    % sign convention: top face = negative of (fn, ft), bottom = positive
    Tn_top_tsl(iS) = -fn;
    Tt_top_tsl(iS) = -ft;
    Tn_bot_tsl(iS) =  fn;
    Tt_bot_tsl(iS) =  ft;
end

% enforce traction-free tip in plotted TSL tractions (consistent with solver)
Tn_top_tsl(end) = 0;  Tt_top_tsl(end) = 0;
Tn_bot_tsl(end) = 0;  Tt_bot_tsl(end) = 0;

% Sort by local abscissa (from CZ mouth)
[xs, idx] = sort(xlocal);

Dn = Dn(idx); Dt = Dt(idx);

Tn_top_bulk = Tn_top_bulk(idx);  Tt_top_bulk = Tt_top_bulk(idx);
Tn_bot_bulk = Tn_bot_bulk(idx);  Tt_bot_bulk = Tt_bot_bulk(idx);

Tn_top_tsl  = Tn_top_tsl(idx);   Tt_top_tsl  = Tt_top_tsl(idx);
Tn_bot_tsl  = Tn_bot_tsl(idx);   Tt_bot_tsl  = Tt_bot_tsl(idx);

% De-duplicate repeated stations (shared vertices) after sorting
tolx = 1e-12 * max(1, da_est);
keep = true(size(xs));
keep(2:end) = abs(diff(xs)) > tolx;

xs = xs(keep); Dn = Dn(keep); Dt = Dt(keep);
Tn_top_bulk = Tn_top_bulk(keep); Tt_top_bulk = Tt_top_bulk(keep);
Tn_bot_bulk = Tn_bot_bulk(keep); Tt_bot_bulk = Tt_bot_bulk(keep);
Tn_top_tsl  = Tn_top_tsl(keep);  Tt_top_tsl  = Tt_top_tsl(keep);
Tn_bot_tsl  = Tn_bot_tsl(keep);  Tt_bot_tsl  = Tt_bot_tsl(keep);

% ------------------ Plots ------------------
figure(5); clf(5)

subplot(3,1,1)
plot(xs, Dn, 'r', 'LineWidth',1); hold on
plot(xs, Dt, 'b', 'LineWidth',1);
grid on
legend('Normal $\bar{d}_n$','Tangential $\bar{d}_t$','Interpreter','latex');
xlabel('$x$ from CZ mouth (local), m','Interpreter','latex');
ylabel('relative opening');
xlim([xs(1), xs(end)]);

subplot(3,1,2); hold on; grid on
plot(xs, Tn_top_bulk, 'r-', 'LineWidth',1);
plot(xs, Tt_top_bulk, 'b-', 'LineWidth',1);
plot(xs, Tn_top_tsl,  'k--','LineWidth',1);
plot(xs, Tt_top_tsl,  'k:','LineWidth',1);
legend('$\bar{\sigma}_n^{\rm bulk}$','$\bar{\sigma}_t^{\rm bulk}$', ...
       '$\bar{\sigma}_n^{\rm TSL}$','$\bar{\sigma}_t^{\rm TSL}$', ...
       'Location','best','Interpreter','latex');
xlabel('$x$ from CZ mouth (local), m','Interpreter','latex');
ylabel({'relative','cohesive traction','(top face)'});
xlim([xs(1), xs(end)]);

subplot(3,1,3); hold on; grid on
plot(xs, Tn_bot_bulk, 'r-', 'LineWidth',1);
plot(xs, Tt_bot_bulk, 'b-', 'LineWidth',1);
plot(xs, Tn_bot_tsl,  'k--','LineWidth',1);
plot(xs, Tt_bot_tsl,  'k:','LineWidth',1);
legend('$\bar{\sigma}_n^{\rm bulk}$','$\bar{\sigma}_t^{\rm bulk}$', ...
       '$\bar{\sigma}_n^{\rm TSL}$','$\bar{\sigma}_t^{\rm TSL}$', ...
       'Location','best','Interpreter','latex');
xlabel('$x$ from CZ mouth (local), m','Interpreter','latex');
ylabel({'relative','cohesive traction','(bottom face)'});
xlim([xs(1), xs(end)]);
end
