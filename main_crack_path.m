function Path = main_crack_path(varargin)
%MAIN_CRACK_PATH  Polyline crack trajectory via repeated calls to main_sweep_theta2.
%
% Concept:
%   Step k:
%     - given physical polyline Pphys^(k), call main_sweep_theta2(Pphys^(k))
%       -> theta_k^*, ell_adv_k
%     - advance the physical tip by ell_adv_k along theta_k^*
%     - append new vertex to Pphys^(k+1)
%     - repeat
%
% Stage-1 version:
%   - Works with the adaptive predictor/bracket/refinement logic
%     in the new main_sweep_theta2.
%   - Passes angle history from previous crack steps.
%   - Keeps stress output.
%   - No old coarse/fine sweep controls.

p = inputParser;

p.addParameter('maxSteps', 60, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('minAdvance', 1e-6, @(x)isnumeric(x)&&isscalar(x)&&x>=0);

% forwarded to main_sweep_theta2
p.addParameter('solverDisplay', 'off', @(s)ischar(s)||isstring(s));
p.addParameter('useWarmStart', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('wantStarSol', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('verbose', true, @(b)islogical(b)&&isscalar(b));

% adaptive angle-search settings
p.addParameter('theta_default_deg', 0, @(x)isnumeric(x)&&isscalar(x));

p.addParameter('init_left_deg',  1.0, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('init_right_deg', 1.0, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('trendBias', 1.0, @(x)isnumeric(x)&&isscalar(x));

p.addParameter('theta_min_deg', -50, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('theta_max_deg',  50, @(x)isnumeric(x)&&isscalar(x));

p.addParameter('expand_step_deg', 0.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('expand_growth',   2.0, @(x)isnumeric(x)&&isscalar(x)&&x>1);
p.addParameter('max_expand_iter', 8,   @(x)isnumeric(x)&&isscalar(x)&&x>=0);

p.addParameter('theta_tol_deg',   0.01, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('dt_tol',          1e-9, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('max_refine_iter', 8,    @(x)isnumeric(x)&&isscalar(x)&&x>=0);

% stress-output controls
p.addParameter('doStressEachStep', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('stressOutDir', 'stress_steps', @(s)ischar(s)||isstring(s));
p.addParameter('nstress', 2, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('stressFact', 10, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('stressDPI', 600, @(x)isnumeric(x)&&isscalar(x));

% adaptive delta controls
p.addParameter('adaptiveDelta', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('deltaFactor', 1.4, @(x)isnumeric(x)&&isscalar(x)&&x>1);
p.addParameter('deltaRelax', 0.5, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('deltaMin', [], @(x)isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
p.addParameter('deltaMax', [], @(x)isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));

p.parse(varargin{:});
S = p.Results;

doStressEachStep = S.doStressEachStep;
stressOutDir     = char(S.stressOutDir);
nstress          = S.nstress;
stressFact       = S.stressFact;
stressDPI        = S.stressDPI;

if doStressEachStep && ~exist(stressOutDir,'dir')
    mkdir(stressOutDir);
end

% initial physical crack from baseline cfg: take all but last point
% (last leg in cfg_min_geom is the cohesive probe)
C0 = cfg_min_geom();
Pphys = C0.Pmid(1:end-1,:);

% baseline probe length
delta0 = norm(C0.Pmid(end,:) - C0.Pmid(end-1,:));

% adaptive-delta bounds
deltaMin = C0.deltaMinFactor * delta0;
deltaMax = C0.deltaMaxFactor * delta0;

% current working trial-leg length
delta = min(deltaMax, max(deltaMin, delta0));

% history containers
Path = struct();
Path.Phist = cell(S.maxSteps+1,1);

Path.theta_star_deg = nan(S.maxSteps,1);
Path.theta_star     = nan(S.maxSteps,1);
Path.sig_star       = nan(S.maxSteps,1);
Path.dt_star        = nan(S.maxSteps,1);
Path.ell_adv        = nan(S.maxSteps,1);
Path.j_star         = nan(S.maxSteps,1);

Path.delta_used   = nan(S.maxSteps,1);
Path.active_ratio = nan(S.maxSteps,1);

Path.bracket   = cell(S.maxSteps,1);
Path.predictor = cell(S.maxSteps,1);
Path.history   = cell(S.maxSteps,1);

% additional diagnostics
Path.delta_used        = nan(S.maxSteps,1);   % actual probe length used at step k
Path.active_ratio      = nan(S.maxSteps,1);   % ell_adv / delta_used

Path.theta_pred_deg    = nan(S.maxSteps,1);   % predictor angle (deg)
Path.theta_pred_err_deg= nan(S.maxSteps,1);   % predictor error (deg)

Path.nExpand           = nan(S.maxSteps,1);   % number of bracket expansions
Path.nRefine           = nan(S.maxSteps,1);   % number of refinement iterations

Path.thetaL_deg        = nan(S.maxSteps,1);   % left bracket boundary
Path.thetaR_deg        = nan(S.maxSteps,1);   % right bracket boundary

Path.tip_xy            = nan(S.maxSteps,2);   % crack tip coordinates
Path.nThetaEval        = nan(S.maxSteps,1);   % total number of theta evaluations at step k

Path.Phist{1} = Pphys;

for k = 1:S.maxSteps

    % angle history for predictor
    theta_prev_deg = [];
    theta_prevprev_deg = [];

    if k >= 2 && isfinite(Path.theta_star_deg(k-1))
        theta_prev_deg = Path.theta_star_deg(k-1);
    end
    if k >= 3 && isfinite(Path.theta_star_deg(k-2))
        theta_prevprev_deg = Path.theta_star_deg(k-2);
    end

    step = main_sweep_theta2(Pphys, ...
        'delta', delta, ...
        'theta_prev_deg', theta_prev_deg, ...
        'theta_prevprev_deg', theta_prevprev_deg, ...
        'theta_default_deg', S.theta_default_deg, ...
        'init_left_deg', S.init_left_deg, ...
        'init_right_deg', S.init_right_deg, ...
        'trendBias', S.trendBias, ...
        'theta_min_deg', S.theta_min_deg, ...
        'theta_max_deg', S.theta_max_deg, ...
        'expand_step_deg', S.expand_step_deg, ...
        'expand_growth', S.expand_growth, ...
        'max_expand_iter', S.max_expand_iter, ...
        'theta_tol_deg', S.theta_tol_deg, ...
        'dt_tol', S.dt_tol, ...
        'max_refine_iter', S.max_refine_iter, ...
        'solverDisplay', S.solverDisplay, ...
        'useWarmStart', S.useWarmStart, ...
        'wantStarSol', S.wantStarSol, ...
        'verbose', S.verbose, ...
        'stepIndex', k);

       % store step outputs
    Path.theta_star_deg(k) = step.theta_star_deg;
    Path.theta_star(k)     = step.theta_star;
    Path.sig_star(k)       = step.sig_star;
    Path.dt_star(k)        = step.dt_star;
    Path.ell_adv(k)        = step.ell_adv_star;
    Path.j_star(k)         = step.j_star;

    Path.delta_used(k) = delta;
    if isfinite(step.ell_adv_star) && isfinite(delta) && delta > 0
        Path.active_ratio(k) = step.ell_adv_star / delta;
    end

    Path.bracket{k}   = step.bracket;
    Path.predictor{k} = step.predictor;
    Path.history{k}   = step.history;

    % ============================================================
    % === Additional diagnostics extraction
    % ============================================================

    % --- delta actually used
    if isfield(step,'delta_used')
        Path.delta_used(k) = step.delta_used;
    else
        Path.delta_used(k) = delta; % fallback
    end

    % --- active cohesive ratio
    if isfinite(step.ell_adv_star) && Path.delta_used(k) > 0
        Path.active_ratio(k) = step.ell_adv_star / Path.delta_used(k);
    end

    % --- predictor info
    if isfield(step,'predictor') && isstruct(step.predictor)
        if isfield(step.predictor,'theta_pred_deg')
            Path.theta_pred_deg(k) = step.predictor.theta_pred_deg;
            Path.theta_pred_err_deg(k) = ...
                step.theta_star_deg - step.predictor.theta_pred_deg;
        end
    end

    % --- bracket info
    if isfield(step,'bracket') && isstruct(step.bracket)
        if isfield(step.bracket,'nExpand')
            Path.nExpand(k) = step.bracket.nExpand;
        end
        if isfield(step.bracket,'thetaL_deg')
            Path.thetaL_deg(k) = step.bracket.thetaL_deg;
        end
        if isfield(step.bracket,'thetaR_deg')
            Path.thetaR_deg(k) = step.bracket.thetaR_deg;
        end
       
        % --- total number of theta evaluations
        if isfield(step,'history') && isstruct(step.history) && isfield(step.history,'theta_deg')
            Path.nThetaEval(k) = numel(step.history.theta_deg);
        end
    end

    % --- refinement info
    if isfield(step,'history') && isstruct(step.history)
        if isfield(step.history,'nRefine')
            Path.nRefine(k) = step.history.nRefine;
        end
    end

    % ============================================================
    % === Stress field output at the selected theta* (critical state)
    % ============================================================
    if doStressEachStep && isfield(step,'sol_star') && isstruct(step.sol_star) && step.sol_star.flag > 0

        hStep = fig_step_summary(step, Path, k, ...
        'nstress', nstress, ...
        'stressFact', stressFact, ...
        'figPosition', [100 100 1500 500]);

        sumname_png = fullfile(stressOutDir, sprintf('summary_step_%03d.png', k));
        sumname_fig = fullfile(stressOutDir, sprintf('summary_step_%03d.fig', k));

        if exist('exportgraphics','file') == 2
            exportgraphics(hStep, sumname_png, 'Resolution', stressDPI);
        else
            print(hStep, sumname_png, '-dpng', sprintf('-r%d', stressDPI));
        end

        savefig(hStep, sumname_fig);
        %close(hStep);

        % Standard stress plotter
        %{

        sol = step.sol_star;

        DispStressExt(sol.mesh, sol.u, sol.cz, sol.mat, sol.quad, nstress, stressFact);

        % Optional title for saved frames
        figure(2);
        title(sprintf('$k=%d,\\ \\theta^*=%.4f^\\circ,\\ \\sigma_{cr}=%.6g$', ...
            k, step.theta_star_deg, step.sig_star), 'Interpreter','latex');

        
        % Save figure
        fname = fullfile(stressOutDir, sprintf('stress_step_%03d.png', k));

        if exist('exportgraphics','file') == 2
            exportgraphics(figure(2), fname, 'Resolution', stressDPI);
        else
            print(figure(2), fname, '-dpng', sprintf('-r%d', stressDPI));
        end

        % --- Save MATLAB figure (editable)
        figfname = fullfile(stressOutDir, sprintf('stress_step_%03d.fig', k));
        savefig(figure(2), figfname);
        %}
    end

    if S.verbose
        fprintf(['k=%d | theta*=%.3f deg | pred=%.3f | err=%.3f | ' ...
            'ratio=%.3f | nThetaEval=%d\n'], ...
            k, ...
            Path.theta_star_deg(k), ...
            Path.theta_pred_deg(k), ...
            Path.theta_pred_err_deg(k), ...
            Path.active_ratio(k), ...
            Path.nThetaEval(k));
    end

    if ~(isfinite(step.ell_adv_star)) || step.ell_adv_star < S.minAdvance
        if S.verbose
            fprintf('Stop: advance too small or invalid.\n');
        end
        break;
    end

    % update physical polyline
    tip = Pphys(end,:);
    th  = step.theta_star; % rad
    Pnew = tip + step.ell_adv_star*[cos(th), sin(th)];
    Pphys = [Pphys; Pnew];

    % store tip position
    Path.tip_xy(k,:) = Pnew;

    Path.Phist{k+1} = Pphys;

    % update delta for the next step
    if S.adaptiveDelta
        delta_target = C0.deltaFactor * step.ell_adv_star;
        delta = (1 - C0.deltaRelax) * delta + C0.deltaRelax * delta_target;
        delta = min(deltaMax, max(deltaMin, delta));
    end
end

% trim unused tail
lastP = find(~cellfun(@isempty, Path.Phist), 1, 'last');
Path.Phist = Path.Phist(1:lastP);

lastK = find(isfinite(Path.theta_star_deg), 1, 'last');
if isempty(lastK)
    lastK = 0;
end

Path.theta_star_deg = Path.theta_star_deg(1:lastK);
Path.theta_star     = Path.theta_star(1:lastK);
Path.sig_star       = Path.sig_star(1:lastK);
Path.dt_star        = Path.dt_star(1:lastK);
Path.ell_adv        = Path.ell_adv(1:lastK);
Path.j_star         = Path.j_star(1:lastK);

Path.delta_used   = Path.delta_used(1:lastK);
Path.active_ratio = Path.active_ratio(1:lastK);

Path.bracket   = Path.bracket(1:lastK);
Path.predictor = Path.predictor(1:lastK);
Path.history   = Path.history(1:lastK);

Path.nThetaEval = Path.nThetaEval(1:lastK);


C0 = cfg_min_geom();
F = make_paper_figures(Path, 'C0', C0, 'outDir', 'paper_figs');

end