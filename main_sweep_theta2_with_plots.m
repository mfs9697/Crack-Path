%% main_sweep_theta2.m
% Unified coarse + fine sweep over theta2 for the critical-state solve (f1_crit3),
% including active cohesive process-zone length extraction.
%
% What it does:
%   (1) COARSE sweep over theta2
%   (2) Local quadratic fit of d_t^{mouth}(theta2) on 3-point neighborhood around min |d_t|
%       -> predicted optimum = root of the fitted quadratic (closest root)
%   (3) FINE sweep in a window around the predicted optimum
%   (4) Final theta* estimate = robust local linear root of d_t^{mouth}(theta2) on fine data
%
% What it plots (per your request):
%   - d_t^{mouth}(theta2) for coarse (with quadratic overlay) + fine
%   - sigma_cr(theta2) for coarse + fine
%   - active process-zone length ell_pz(theta2) for coarse + fine
%
% Saves:
%   coarse_sweep_theta2.mat, fine_sweep_theta2.mat (struct fields)

clear; clc;

%% ---------------- USER SETTINGS ----------------
% Coarse sweep set
Theta2_deg_coarse = -5:1:2;

% Fine sweep construction
fine_half_width_deg = .7;     % half-width around predicted optimum (deg)
dtheta_fine_deg     = 0.1;    % fine step (deg)

% Warm start & solver verbosity
useWarmStart   = true;
solverDisplay  = 'off';        % 'off' | 'iter' | 'iter-detailed'

% Save
saveCoarse = 'coarse_sweep_theta2.mat';
saveFine   = 'fine_sweep_theta2.mat';

%% ==================== COARSE SWEEP ====================
coarse = run_theta2_sweep(Theta2_deg_coarse, solverDisplay);
save(saveCoarse, '-struct', 'coarse');
fprintf('Saved: %s\n', saveCoarse);

% Local quadratic fit on dt(theta) near min |dt| and predicted root (closest)
fitC = local_quadratic_fit(coarse.Theta2_deg, coarse.dt_mouth, coarse.flag);
theta_pred_deg = fitC.theta_root_deg;

fprintf('\n==== COARSE FIT ====\n');
fprintf('theta_pred (root of local quadratic dt fit) = %.6f deg\n', theta_pred_deg);

% Plot coarse dt with quadratic overlay (only dt plot needs parabola)
plot_dt_with_quadratic(coarse, fitC, 'Coarse sweep: d_t^{mouth}(\theta_2) with local quadratic fit');

% Plot coarse sigma and process-zone length
plot_sweep_basic(coarse, 'Coarse sweep');

%% ==================== BUILD FINE SET ====================
thetaL = theta_pred_deg - fine_half_width_deg;
thetaR = theta_pred_deg + fine_half_width_deg;

Theta2_deg_fine = thetaL:dtheta_fine_deg:thetaR;
if Theta2_deg_fine(end) < thetaR
    Theta2_deg_fine(end+1) = thetaR;
end

fprintf('\n==== FINE SWEEP SETUP ====\n');
fprintf('Window: [%.3f, %.3f] deg, step = %.3f deg, N = %d\n', ...
    thetaL, thetaR, dtheta_fine_deg, numel(Theta2_deg_fine));

%% ==================== FINE SWEEP ====================
fine = run_theta2_sweep(Theta2_deg_fine, solverDisplay);
save(saveFine, '-struct', 'fine');
fprintf('Saved: %s\n', saveFine);

% Final optimum estimate (recommended): root of dt(theta) from fine data
theta_star_deg = estimate_root_local_linear(fine.Theta2_deg, fine.dt_mouth, fine.flag);

fprintf('\n==== FINAL THETA* (fine sweep) ====\n');
fprintf('theta* (dt=0 root) = %.6f deg\n', theta_star_deg);

% Plot fine sweep results (dt, sigma, process-zone length)
plot_sweep_basic(fine, 'Fine sweep');

%% ========================= FUNCTIONS =========================

function S = run_theta2_sweep(Theta2_deg, solverDisplay)
%RUN_THETA2_SWEEP Unified sweep routine for both coarse and fine sweeps.
%
% Returns struct S with fields:
%   Theta2_deg, Theta2, sig_cr, dt_mouth, ell_pz, flag, iters, Fnorm

    Theta2 = Theta2_deg*pi/180;
    nT = numel(Theta2);

    sig_cr   = nan(nT,1);
    dt_mouth = nan(nT,1);
    ell_pz   = nan(nT,1);   % active process-zone length
    flag     = nan(nT,1);
    iters    = nan(nT,1);
    Fnorm    = nan(nT,1);

    for i = 1:nT
        theta2 = Theta2(i);
        fprintf('\n==== SWEEP: i=%d/%d, theta2 = %+g deg ====\n', i, nT, Theta2_deg(i));

        %% config + geometry
        C = cfg_min_geom(theta2);
        [mesh, ~, ~, elod, Stif, cohes] = geom_pencil(C);

        % Cohesive-zone bundle passed to f1_crit3 and tsl2
        cz = struct('cohes',cohes, 'Q',C.Q,...
            'Dmax',C.Dmax, 'sigmax',C.sigmax, ...
            'a1',C.a1, 'a2',C.a2, 'rphi',C.rphi);

        ndof = 2*size(mesh.coord,1);

        %% Jacobian sparsity pattern
        Jpat = spones(Stif);

        mcoh      = size(cohes,1);
        ncoh_elem = (mcoh-1)/2;

        for q = 1:ncoh_elem
            ind1 = 2*(q-1)+(1:3);
            up   = cohes(ind1,1);
            lo   = cohes(ind1,2);
            dof  = reshape([2*up-1; 2*up; 2*lo-1; 2*lo], [], 1);
            Jpat(dof,dof) = 1;
        end

        Jpat_aug = blkdiag(Jpat, speye(1));

        %% fsolve options
        x_scale = max(1e-9, (C.sig_guess/max(cz.sigmax(1),1))*C.B) * ones(ndof+1,1);

        opts_solve = optimoptions('fsolve', ...
            'Algorithm','trust-region-dogleg', ...
            'SpecifyObjectiveGradient', true, ...
            'JacobPattern', Jpat_aug, ...
            'ScaleProblem','jacobian', ...
            'TypicalX', x_scale, ...
            'FunctionTolerance',1e-12, ...
            'StepTolerance',1e-12, ...
            'OptimalityTolerance',1e-12, ...
            'MaxIterations',100, ...
            'Display',solverDisplay);

        %% initial guess
        x0 = [zeros(ndof,1); C.sig_guess];

        fun = @(xx) f1_crit3(xx, mesh, cz, Stif, elod, cohes);

        [x, fval, flg, out] = fsolve(fun, x0, opts_solve);

        flag(i)  = flg;
        iters(i) = out.iterations;
        Fnorm(i) = norm(fval);

        u_sol     = x(1:ndof);
        sig_cr(i) = x(end);

        %% mouth opening
        i_mouth = 1;
        [~, Dt_mouth] = coh_opening_at_station(u_sol, cohes, cz.Q, i_mouth);
        dt_mouth(i) = Dt_mouth;

        %% active process-zone length (along last leg)
        % delta is cohesive extension length (last polyline leg length)
        delta = C.L(end);

        ell_pz(i) = active_process_zone_length(u_sol, cohes, cz.Q, cz, delta);

        fprintf('flag=%d, iters=%d, ||F||=%.3e, sig=%.6g, Dt_mouth=%.6g, ell_pz/delta=%.6g\n', ...
            flag(i), iters(i), Fnorm(i), sig_cr(i), dt_mouth(i), ell_pz(i));
    end

    S.Theta2_deg = Theta2_deg(:);
    S.Theta2     = Theta2(:);
    S.sig_cr     = sig_cr(:);
    S.dt_mouth   = dt_mouth(:);
    S.ell_pz     = ell_pz(:);
    S.flag       = flag(:);
    S.iters      = iters(:);
    S.Fnorm      = Fnorm(:);
end

function ell_pz = active_process_zone_length(u, cohes, Q, cz, delta)
%ACTIVE_PROCESS_ZONE_LENGTH (opening-based)
% Defines active CZ length using xi_n = d_n/d_nmax compared to a1.
%
% Active zone criterion (tip-attached contiguous):
%   xi_n(s) >= a1  is considered active,
% and ell_pz is the contiguous length near the TIP satisfying this criterion.

    mcoh      = size(cohes,1);
    ncoh_elem = (mcoh-1)/2;

    % station spacing (vertex-mid-vertex...): total stations = 2*ncoh_elem + 1
    ds = delta / (2*ncoh_elem);

    dnmax = cz.Dmax(1);
    a1    = cz.a1;

    xi_n = nan(mcoh,1);

    for i = 1:mcoh
        [Dn, ~] = coh_opening_at_station(u, cohes, Q, i);  % (Dn,Dt) returned
        xi_n(i) = Dn / dnmax;
    end

    % Determine contiguous active zone adjacent to MOUTH:
    % Find first station from the MOUTH where xi_n < a1.
    % Everything after that is inactive.
    j_first_inactive = find(xi_n < a1, 1, 'first');

    ell_pz = max(0, (j_first_inactive-1) * ds/delta);
end

function [fitC] = local_quadratic_fit(theta_deg, dt, flag)
%LOCAL_QUADRATIC_FIT fit dt(theta)=a theta^2 + b theta + c near min |dt|.
    ok = (flag>0) & isfinite(dt);
    th = theta_deg(ok);
    y  = dt(ok);

    [~, j0] = min(abs(y));
    if j0==1 || j0==numel(y)
        error('local_quadratic_fit: minimizer at boundary; widen sweep.');
    end

    x3 = th([j0-1 j0 j0+1]);
    y3 = y ([j0-1 j0 j0+1]);

    p = polyfit(x3, y3, 2); % [a b c]
    a = p(1); b = p(2); c = p(3);

    r = roots(p);
    x0 = th(j0);
    [~, ir] = min(abs(r - x0));
    theta_root_deg = real(r(ir));

    fitC.a = a; fitC.b = b; fitC.c = c;
    fitC.x3 = x3; fitC.y3 = y3;
    fitC.theta_root_deg = theta_root_deg;
end

function plot_dt_with_quadratic(S, fitC, plotTitle)
%PLOT_DT_WITH_QUADRATIC plots dt_mouth(theta) and overlays local quadratic fit.
    ok = (S.flag>0) & isfinite(S.dt_mouth);
    th = S.Theta2_deg(ok);
    y  = S.dt_mouth(ok);
    [th, I] = sort(th); y = y(I);

    figure; hold on; box on;
    plot(th, y, 'o-','LineWidth',1.5);
    yline(0,'-','LineWidth',1.0);

    xx = linspace(min(fitC.x3)-2, max(fitC.x3)+2, 200);
    yy = fitC.a*xx.^2 + fitC.b*xx + fitC.c;
    plot(xx, yy, 'LineWidth',1.5);

    xline(fitC.theta_root_deg,'--','LineWidth',1.2);

    xlabel('$\theta_2$ [deg]','Interpreter','latex');
    ylabel('$d_t^{\mathrm{mouth}}$','Interpreter','latex');
    title(plotTitle);
    legend('data','0-line','local quadratic','predicted root','Location','best');
    grid on;
end

function plot_sweep_basic(S, baseTitle)
%PLOT_SWEEP_BASIC
% Plots:
%   (1) d_t^{mouth}(\theta_2)
%   (2) \sigma_{cr}(\theta_2)
%   (3) \ell_{pz}(\theta_2)
% using LaTeX interpreter for all text.

    ok = (S.flag>0) & isfinite(S.dt_mouth) & isfinite(S.sig_cr) & isfinite(S.ell_pz);
    th = S.Theta2_deg(ok);
    dt = S.dt_mouth(ok);
    sg = S.sig_cr(ok);
    lp = S.ell_pz(ok);

    [th, I] = sort(th);
    dt = dt(I);
    sg = sg(I);
    lp = lp(I);

    % ---- 1) d_t^{mouth}(\theta_2) ----
    figure; hold on; box on;
    plot(th, dt, 'o-','LineWidth',1.5);
    yline(0,'-','LineWidth',1.0);

    xlabel('$\theta_2\,[\mathrm{deg}]$','Interpreter','latex');
    ylabel('$d_t^{\mathrm{mouth}}$','Interpreter','latex');
    title([baseTitle, ': $d_t^{\mathrm{mouth}}(\theta_2)$'], ...
        'Interpreter','latex');

    grid on;

    % ---- 2) \sigma_{cr}(\theta_2) ----
    figure; hold on; box on;
    plot(th, sg, 'o-','LineWidth',1.5);

    xlabel('$\theta_2\,[\mathrm{deg}]$','Interpreter','latex');
    ylabel('$\sigma_{\mathrm{cr}}$','Interpreter','latex');
    title([baseTitle, ': $\sigma_{\mathrm{cr}}(\theta_2)$'], ...
        'Interpreter','latex');

    grid on;

    % ---- 3) \ell_{pz}(\theta_2) ----
    figure; hold on; box on;
    plot(th, lp, 'o-','LineWidth',1.5);

    xlabel('$\theta_2\,[\mathrm{deg}]$','Interpreter','latex');
    ylabel('$\ell_{\mathrm{pz}}$','Interpreter','latex');
    title([baseTitle, ': active process-zone length $\ell_{\mathrm{pz}}(\theta_2)$'], ...
        'Interpreter','latex');

    grid on;
end


function theta_root_deg = estimate_root_local_linear(theta_deg, dt, flag)
%ESTIMATE_ROOT_LOCAL_LINEAR robust dt=0 root using nearest sign change.
    ok = (flag>0) & isfinite(dt);
    th = theta_deg(ok);
    y  = dt(ok);

    [th, I] = sort(th);
    y = y(I);

    [~, j0] = min(abs(y));

    jL = []; jR = [];
    for j = j0:-1:2
        if y(j-1)*y(j) <= 0
            jL = j-1; jR = j; break;
        end
    end
    if isempty(jL)
        for j = j0:1:numel(y)-1
            if y(j)*y(j+1) <= 0
                jL = j; jR = j+1; break;
            end
        end
    end
    if isempty(jL)
        error('estimate_root_local_linear: no sign change found.');
    end

    x1 = th(jL); x2 = th(jR);
    y1 = y(jL);  y2 = y(jR);

    if abs(y2-y1) < eps
        theta_root_deg = 0.5*(x1+x2);
    else
        theta_root_deg = x1 - y1*(x2-x1)/(y2-y1);
    end
end

function [Dn, Dt] = coh_opening_at_station(u, cohes, Q, i_station)
%COH_OPENING_AT_STATION
% Convention (per your note):
%   D_local = Q * (u_up - u_lo)
%   Dt = D_local(1)
%   Dn = D_local(2)
    up = cohes(i_station,1);
    lo = cohes(i_station,2);

    du = [ u(2*up-1) - u(2*lo-1);
           u(2*up  ) - u(2*lo  ) ];

    if ndims(Q)==3
        Qi = Q(:,:,i_station);
    else
        Qi = Q;
    end

    Dloc = Qi * du;

    Dt = Dloc(1);
    Dn = Dloc(2);
end
