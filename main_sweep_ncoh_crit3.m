function out = main_sweep_ncoh_crit3(varargin)
% MAIN_SWEEP_NCOH_CRIT3  Interface-refinement convergence study for
% sigma_cr vs n_coh using the current cfg() + geom_pencil() + f1_crit3().
%
% Goal:
%   Produce a clean, reproducible convergence table for a representative
%   kink configuration (fixed theta_2), demonstrating robustness with
%   respect to cohesive interface discretization n_coh.
%
% Current solver path:
%   C = cfg();
%   [V, mesh, mat, quad, elod, Stif, cohes] = geom_pencil(C);
%   [F_aug, J_aug] = f1_crit3(x, mesh, cz, Stif, elod, cohes);
%
% Unknown vector:
%   x = [u; sig], where u stacks 2 DOFs per node (ux, uy).
%
% Reported energetic characteristic:
%   hatPsi - c_n
%
% Example:
%   out = main_sweep_ncoh_crit3( ...
%       'theta2_deg', 25, ...
%       'ncoh', [10 20 40 80], ...
%       'hgrad', 1.10, ...
%       'hmax_ratio', 80, ...
%       'save_dir', 'results_ncoh_crit3');

% ---------------------- parameters ----------------------
p = inputParser; p.KeepUnmatched = true;

p.addParameter('ncoh',        [20 40 80],  @(v)isnumeric(v)&&isvector(v));
p.addParameter('theta2_deg',  20,             @(x)isnumeric(x)&&isscalar(x));
p.addParameter('hgrad',       1.10,           @(x)isnumeric(x)&&isscalar(x));
p.addParameter('hmax_ratio',  80,             @(x)isnumeric(x)&&isscalar(x));
p.addParameter('da_ratio',    .06,             @(x)isnumeric(x)&&isscalar(x));
p.addParameter('save_dir',    'results_ncoh_crit3', @(s)ischar(s)||isstring(s));
p.addParameter('sig_guess',   [],             @(x)isnumeric(x)&&isscalar(x));
p.addParameter('display',     'none',         @(s)ischar(s)||isstring(s));
p.addParameter('maxIter',     200,            @(x)isnumeric(x)&&isscalar(x));
p.addParameter('maxFunEvals', 400,            @(x)isnumeric(x)&&isscalar(x));

p.parse(varargin{:});
P = p.Results;
P.save_dir = char(P.save_dir);
if ~exist(P.save_dir,'dir'); mkdir(P.save_dir); end

% ---------------------- base config ----------------------
C0 = cfg();

% Representative fixed kink angle
C0.theta_deg = P.theta2_deg;
C0.alf2      = deg2rad(P.theta2_deg);

% Rebuild local frame at leg 2
calf = cos(C0.alf2); salf = sin(C0.alf2);
C0.malf = [calf, salf; -salf, calf];

% Mesh controls
C0.hgrad = P.hgrad;
C0.hmax  = C0.B / P.hmax_ratio;

% Optional cohesive probe length override
if ~isempty(P.da_ratio)
    C0.da = P.da_ratio * C0.a;
end

% Rebuild geometric points after any angle/da change
C0.V0 = [0, 0];
C0.V1 = C0.a  * [cos(C0.alf1), sin(C0.alf1)];
C0.V2 = C0.V1 + C0.da * [cos(C0.alf2), sin(C0.alf2)];
C0.Pmid = [C0.V0; C0.V1; C0.V2];
C0.holes = {};
C0.hmax_ratio = P.hmax_ratio;

% Optional sigma guess override
if ~isempty(P.sig_guess)
    C0.sig_guess = P.sig_guess;
end

% Disable plots for sweep runs
%C0.plotMesh = 0;

nlist = sort(P.ncoh(:).');
nlev  = numel(nlist);

rec = struct( ...
    'ncoh',[], 'hmin',[], 'dofs',[], 'sigma_cr',[], ...
    'Dn_mouth',[], 'Dt_mouth',[], ...
    'psi_hat',[], 'Psi_hat_minus_cn',[], ...
    'iters',[], 'funcCount',[], 'exitflag',[], 'msg',[], ...
    't_mesh',[], 't_prep',[], 't_solve',[], 't_total',[], ...
    'ratio_mesh_solve',[], 'ratio_prep_solve',[]);
rec(nlev).ncoh = [];

fprintf('\n Interface-refinement sweep with f1_crit3\n');
fprintf(' theta_2 = %.1f deg, h_grad = %.2f, B/h_max = %.1f\n', ...
    P.theta2_deg, P.hgrad, P.hmax_ratio);
fprintf('%8s %10s %8s %14s %14s %s\n', ...
    'n_coh','DOFs','iters','sigma_cr','hatPsi-c_n','exitflag');
fprintf('%s\n', repmat('-',1,88));

% ---------------------- loop over n_coh ----------------------
for k = 1:nlev
    C = C0;
    C.ncoh = nlist(k);
    C.chw  = (1/8) * C.da / C.ncoh;
    hmin_k = C.da / C.ncoh;

    try
        % ---------- (A) Mesh + stiffness ----------
        t0 = tic;
        [mesh, ~, ~, elod, Stif, cohes] = geom_pencil(C);
        t_mesh = toc(t0);

        % ---------- (B) Prep ----------
        t1 = tic;

        cz = build_cz_from_geom_and_cfg(cohes, C);

        ndof = size(mesh.coord,1)*2;
        x0   = zeros(ndof+1,1);
        if isfield(C,'sig_guess') && ~isempty(C.sig_guess) && (C.sig_guess > 0)
            x0(end) = C.sig_guess;
        else
            x0(end) = 0.1;
        end

        % Jacobian sparsity pattern
        Jpat = spones(Stif);

        mcoh = size(cz.cohes,1);
        ncoh_local = (mcoh-1)/2;
        for q = 1:ncoh_local
            ind1 = 2*(q-1) + (1:3);
            up  = cz.cohes(ind1,1);
            lo  = cz.cohes(ind1,2);
            dof = reshape([2*up-1; 2*up; 2*lo-1; 2*lo], [], 1);
            Jpat(dof,dof) = 1;
        end
        Jpat_aug = blkdiag(Jpat, speye(1));

        sig_guess = x0(end);
        sig_ref = max(C.sigmax(1), 1);
        x_scale = max(1e-9, (sig_guess/sig_ref)*C.B) * ones(ndof+1,1);

        opts = optimoptions('fsolve', ...
            'Algorithm','trust-region-dogleg', ...
            'SpecifyObjectiveGradient', true, ...
            'JacobPattern', Jpat_aug, ...
            'ScaleProblem','jacobian', ...
            'TypicalX', x_scale, ...
            'FunctionTolerance',1e-12, ...
            'StepTolerance',1e-12, ...
            'OptimalityTolerance',1e-12, ...
            'MaxIterations', P.maxIter, ...
            'MaxFunctionEvaluations', P.maxFunEvals, ...
            'Display', char(P.display));

        fun = @(x) f1_crit3(x, mesh, cz, Stif, elod, cohes);

        t_prep = toc(t1);

        % ---------- (C) Nonlinear solve ----------
        t2 = tic;
        [x, ~, exitflag, outfs] = fsolve(fun, x0, opts);
        t_solve = toc(t2);

        t_total = t_mesh + t_prep + t_solve;
        sigma_cr = x(end);
        u = x(1:ndof);

        % ---------- (D) Mouth energetic characteristic ----------
        iu = cz.cohes(1,1);
        il = cz.cohes(1,2);

        Uu = [u(2*iu-1); u(2*iu)];
        Ul = [u(2*il-1); u(2*il)];

        Dl   = cz.Q * (Uu - Ul);   % [Dt; Dn]
        Dt_m = Dl(1);
        Dn_m = Dl(2);

        if abs(Dn_m - cz.Dmax(1)) > 1e-6
            error('Critical constraint violated: Dn_mouth != Dn_max (%.3e)', ...
                abs(Dn_m - cz.Dmax(1)));
        end

        xi_t = abs(Dt_m) / cz.Dmax(2);

        cn = (3 - 2*cz.a1 + 3*cz.a2(1))/6;
        phi_n = cz.phi(1);
        phi_t = cz.phi(2);

        [F2, ~, ~] = g_base(xi_t, cz.a1, cz.a2(2), 0);

        psi_hat = cn + (phi_t/phi_n - cn*cz.rphi*sqrt(phi_t/phi_n)) * F2;
        Psi_hat_minus_cn = psi_hat - cn;

        % ---------- bookkeeping ----------
        rec(k).ncoh             = C.ncoh;
        rec(k).hmin             = hmin_k;
        rec(k).dofs             = ndof;
        rec(k).sigma_cr         = sigma_cr;
        rec(k).Dn_mouth         = Dn_m;
        rec(k).Dt_mouth         = Dt_m;
        rec(k).psi_hat          = psi_hat;
        rec(k).Psi_hat_minus_cn = Psi_hat_minus_cn;
        rec(k).iters            = outfs.iterations;
        rec(k).funcCount        = outfs.funcCount;
        rec(k).exitflag         = exitflag;
        rec(k).msg              = 'ok';

        rec(k).t_mesh           = t_mesh;
        rec(k).t_prep           = t_prep;
        rec(k).t_solve          = t_solve;
        rec(k).t_total          = t_total;
        rec(k).ratio_mesh_solve = t_mesh / max(t_solve, eps);
        rec(k).ratio_prep_solve = (t_mesh + t_prep) / max(t_solve, eps);

        fprintf('%8d %10d %8d %14.6f %14.5e  %3d\n', ...
            C.ncoh, ndof, outfs.iterations, sigma_cr, Psi_hat_minus_cn, exitflag);

    catch ME
        rec(k).ncoh             = C.ncoh;
        rec(k).hmin             = hmin_k;
        rec(k).dofs             = NaN;
        rec(k).sigma_cr         = NaN;
        rec(k).Dn_mouth         = NaN;
        rec(k).Dt_mouth         = NaN;
        rec(k).psi_hat          = NaN;
        rec(k).Psi_hat_minus_cn = NaN;
        rec(k).iters            = NaN;
        rec(k).funcCount        = NaN;
        rec(k).exitflag         = NaN;
        rec(k).msg              = ME.message;

        rec(k).t_mesh           = NaN;
        rec(k).t_prep           = NaN;
        rec(k).t_solve          = NaN;
        rec(k).t_total          = NaN;
        rec(k).ratio_mesh_solve = NaN;
        rec(k).ratio_prep_solve = NaN;

        warning('n_coh=%d failed: %s', C.ncoh, ME.message);
        fprintf('%8d %10s %8s %14s %14s  %s\n', ...
            C.ncoh, '—', '—', '—', '—', 'FAILED');
    end

    clear V mesh Stif elod cohes cz
end

fprintf('%s\n', repmat('-',1,88));

% ---------------------- Richardson extrapolation ----------------------
okmask = ~isnan([rec.sigma_cr]) & ([rec.exitflag] > 0);
R      = rec(okmask);

out = struct('records',rec,'h',[],'sigma_cr',[],'p',NaN,'sigma_inf',NaN,'err_rel',NaN);

if numel(R) >= 3
    h   = [R.hmin];
    sig = [R.sigma_cr];

    [h, idx] = sort(h, 'descend');
    sig = sig(idx);

    h3   = h(end-2:end);
    sig3 = sig(end-2:end);

    r1 = h3(1)/h3(2);
    r2 = h3(2)/h3(3);

    if abs(r1 - r2) <= 1e-14
        r  = r1;
        s0 = sig3(1); s1 = sig3(2); s2 = sig3(3);

        num = (s0 - s1);
        den = (s1 - s2);

        if abs(den) > eps && abs(num/den) > 0
            p_est = log(abs(num/den)) / log(r);
            sigma_inf = s2 + (s2 - s1) / (r^p_est - 1);
            err_rel   = abs(s2 - sigma_inf) / max(1, abs(sigma_inf));

            out.p         = p_est;
            out.sigma_inf = sigma_inf;
            out.err_rel   = err_rel;

            fprintf('Richardson: sigma_inf = %.6e, rel.err = %.2e\n', ...
                out.sigma_inf, out.err_rel);
        else
            warning('Degenerate differences; skipping Richardson estimate.');
        end
    else
        warning('Non-uniform refinement ratios; skipping Richardson estimate.');
    end

    out.h        = h;
    out.sigma_cr = sig;
else
    warning('Need at least 3 successful levels for Richardson analysis.');
end

% ---------------------- save artifacts ----------------------
try
    T = struct2table(rec);
    writetable(T, fullfile(P.save_dir,'sweep_ncoh_crit3.csv'));

    figure(101); clf
    loglog(out.h, out.sigma_cr, 'o-','LineWidth',1.2,'MarkerSize',6); grid on
    set(gca,'XDir','reverse');
    xlabel('$h_{\min} = d a / n_{\mathrm{coh}}$','Interpreter','latex');
    ylabel('$\sigma_{\mathrm{cr}}$','Interpreter','latex');
    title(sprintf('Critical load vs. interface resolution with f1\\_crit3 (\\theta_2=%.1f^\\circ)', ...
        P.theta2_deg));

    if ~isnan(out.sigma_inf)
        hold on
        yline(out.sigma_inf,'--','\sigma_\infty','LabelVerticalAlignment','bottom');
        text(out.h(end), out.sigma_cr(end), ...
            sprintf('  rel.err \\approx %.1e', out.err_rel), ...
            'VerticalAlignment','bottom');
    end

    saveas(gcf, fullfile(P.save_dir,'sigma_vs_hmin_crit3.png'));
    save(fullfile(P.save_dir,'sweep_ncoh_crit3.mat'),'rec','out','P');
catch
    % best-effort
end

fprintf('Saved results to %s\n', P.save_dir);

end

% ============================================================
function cz = build_cz_from_geom_and_cfg(cohes, C)
% BUILD_CZ_FROM_GEOM_AND_CFG
% V      : [3x2] = [mouth; knee; tip]
% cohes  : [(2K+1)x2] node pairs (upper, lower), ordered knee..tip
% C      : cfg() struct with Dmax, sigmax, a1, a2, rphi, phi, malf

cz = struct();
cz.cohes = cohes;

% Current project convention
cz.Q = C.malf;
cz.m = C.malf;

% Cohesive parameters
cz.phi    = C.phi(1:2);
cz.Dmax   = C.Dmax(1:2);      % [Dn_max, Dt_max]
cz.sigmax = C.sigmax(1:2);    % [sig_n_max, sig_t_max]
cz.a1     = C.a1;
cz.a2     = C.a2;
cz.rphi   = C.rphi;
end