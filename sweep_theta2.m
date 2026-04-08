function out = sweep_theta2(Pphys, varargin)
%MAIN_SWEEP_THETA2  One growth step: adaptive search for the kinking angle.
%
% Stage-1 version:
%   - Uses solve_one_theta(...) as the single source of truth.
%   - Replaces fixed coarse/fine sweeps by:
%       predictor -> initial window -> adaptive bracket -> bracket refinement
%   - Keeps only minimal command-window output.
%   - No plotting in this version.
%
% Inputs:
%   Pphys : [Np x 2] physical crack polyline (mouth -> current tip).
%
% Name-value options:
%   'cfg'                     : optional configuration override struct
%   'x0_prev'                 : optional warm start for the initial center trial
%
%   'delta'                   : cohesive probe length
%   'theta_prev_deg'          : previously accepted kinking angle (deg)
%   'theta_prevprev_deg'      : angle from two steps ago (deg)
%   'theta_default_deg'       : default predictor if no history (deg, default 0)
%
%   'init_left_deg'           : initial window left half-width  (default 1.0)
%   'init_right_deg'          : initial window right half-width (default 0.5)
%   'trendBias'               : extrapolation factor beta for predictor (default 1.0)
%
%   'theta_min_deg'           : lower angular bound (default -50)
%   'theta_max_deg'           : upper angular bound (default  50)
%
%   'expand_step_deg'         : first expansion step (default 0.5)
%   'expand_growth'           : geometric growth factor for expansion (default 2.0)
%   'max_expand_iter'         : max bracket expansions (default 8)
%
%   'theta_tol_deg'           : refinement tolerance in theta (default 0.02)
%   'dt_tol'                  : refinement tolerance in dt_mouth (default 1e-8)
%   'max_refine_iter'         : max refinement steps (default 8)
%
%   'solverDisplay'           : fsolve display
%   'useWarmStart'            : true/false
%   'wantStarSol'             : true/false
%   'verbose'                 : true/false
%   'stepIndex'               : optional crack-step index for printing
%
% Output:
%   out.theta_star_deg, out.theta_star
%   out.sig_star, out.dt_star
%   out.ell_adv_star, out.j_star
%   out.Pmid_star
%   out.sol_star        (optional)
%   out.history         evaluated trials
%   out.bracket         bracket information
%   out.predictor       predictor information
%   out.cfg_used        configuration actually used

% ---------------- parse inputs ----------------
p = inputParser;
p.addRequired('Pphys', @(x)isnumeric(x)&&size(x,2)==2&&size(x,1)>=2);

p.addParameter('cfg', [], @(x) isempty(x) || isstruct(x));
p.addParameter('x0_prev', [], @(x) isempty(x) || isnumeric(x));

p.addParameter('delta', [], @(x)isnumeric(x)&&isscalar(x)&&x>0);

p.addParameter('theta_prev_deg', [], @(x)isempty(x) || (isnumeric(x)&&isscalar(x)));
p.addParameter('theta_prevprev_deg', [], @(x)isempty(x) || (isnumeric(x)&&isscalar(x)));
p.addParameter('theta_default_deg', 0, @(x)isnumeric(x)&&isscalar(x));

p.addParameter('init_left_deg',  1.0, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('init_right_deg', 1.0, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('trendBias', 1.0, @(x)isnumeric(x)&&isscalar(x));

p.addParameter('theta_min_deg', -50, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('theta_max_deg',  50, @(x)isnumeric(x)&&isscalar(x));

p.addParameter('expand_step_deg', 0.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('expand_growth',   2.0, @(x)isnumeric(x)&&isscalar(x)&&x>1);
p.addParameter('max_expand_iter', 8,   @(x)isnumeric(x)&&isscalar(x)&&x>=0);

p.addParameter('theta_tol_deg',   0.02, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('dt_tol',          1e-8, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('max_refine_iter', 8,    @(x)isnumeric(x)&&isscalar(x)&&x>=0);

p.addParameter('solverDisplay', 'off', @(s)ischar(s)||isstring(s));
p.addParameter('useWarmStart', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('wantStarSol',  true, @(b)islogical(b)&&isscalar(b));
p.addParameter('verbose',      true, @(b)islogical(b)&&isscalar(b));
p.addParameter('stepIndex', [], @(x)isempty(x) || (isnumeric(x)&&isscalar(x)));

p.parse(Pphys, varargin{:});
S = p.Results;

solverDisplay = char(S.solverDisplay);

% ---------------- base cfg + probe length ----------------
if isempty(S.cfg)
    C0 = cfg_min_geom();
else
    C0 = S.cfg;
end

if isempty(S.delta)
    if isfield(C0, 'Pmid') && size(C0.Pmid,1) >= 2
        S.delta = norm(C0.Pmid(end,:) - C0.Pmid(end-1,:));
    elseif isfield(C0, 'delta') && isnumeric(C0.delta) && isscalar(C0.delta) && C0.delta > 0
        S.delta = C0.delta;
    else
        error('main_sweep_theta2:MissingDelta', ...
            ['No probe length provided. Pass ''delta'' explicitly or ', ...
             'provide it in cfg as C0.delta or through C0.Pmid.']);
    end
end
dLeg = S.delta;

% ---------------- predictor ----------------
pred = make_theta_predictor( ...
    S.theta_prev_deg, S.theta_prevprev_deg, ...
    S.theta_default_deg, S.trendBias, ...
    S.theta_min_deg, S.theta_max_deg);

theta_pred_deg = pred.theta_pred_deg;

% Initial asymmetric window; if recent trend is positive, flip asymmetry
if isfinite(pred.trend_deg) && pred.trend_deg > 0
    dL0 = S.init_right_deg;
    dR0 = S.init_left_deg;
else
    dL0 = S.init_left_deg;
    dR0 = S.init_right_deg;
end

thetaL0_deg = max(S.theta_min_deg, theta_pred_deg - dL0);
thetaC0_deg = min(max(theta_pred_deg, S.theta_min_deg), S.theta_max_deg);
thetaR0_deg = min(S.theta_max_deg, theta_pred_deg + dR0);

tip = Pphys(end,:);

if S.verbose
    if ~isempty(S.stepIndex)
        fprintf('\n=================== CRACK STEP %d ===================\n', S.stepIndex);
    else
        fprintf('\n=================== CRACK STEP ===================\n');
    end

    fprintf('tip               = [%g, %g]\n', tip(1), tip(2));
    fprintf('theta_pred        = %8.4f deg\n', theta_pred_deg);
    fprintf('initial window    = [%8.4f, %8.4f, %8.4f] deg\n', ...
        thetaL0_deg, thetaC0_deg, thetaR0_deg);
    fprintf('delta             = %.6g,  delta/a= %.6g \n', dLeg,dLeg/C0.a);

    if isfield(C0,'ncoh') && ~isempty(C0.ncoh)
        fprintf('ncoh              = %d\n', C0.ncoh);
    end
    if isfield(C0,'a1') && ~isempty(C0.a1)
        fprintf('a1                = %.6g\n', C0.a1);
    end
    if isfield(C0,'sigmax') && numel(C0.sigmax) >= 2
        fprintf('sigmax_n, sigmax_t= %.6g, %.6g\n', C0.sigmax(1), C0.sigmax(2));
    end

    %fprintf('note              = beta == ell_adv\n');
end

% ---------------- trial cache ----------------
cache = struct();
cache.theta_deg = [];
cache.trials    = {};

% evaluate initial triplet
% Use explicit x0_prev only for the initial center trial, if supplied.
[cache, iL] = eval_trial_cached(cache, thetaL0_deg, Pphys, C0, dLeg, [],         solverDisplay, S.useWarmStart, S.verbose); %#ok<NASGU>
[cache, iC] = eval_trial_cached(cache, thetaC0_deg, Pphys, C0, dLeg, S.x0_prev,  solverDisplay, S.useWarmStart, S.verbose); %#ok<NASGU>
[cache, iR] = eval_trial_cached(cache, thetaR0_deg, Pphys, C0, dLeg, [],         solverDisplay, S.useWarmStart, S.verbose); %#ok<NASGU>

% ---------------- bracket search ----------------
[bracket, cache] = find_bracket_adaptive(cache, Pphys, C0, dLeg, S, solverDisplay);

% ---------------- refinement / fallback ----------------
if bracket.found
    [trial_star, cache] = refine_bracketed_root(bracket, cache, Pphys, C0, dLeg, S, solverDisplay);
else
    trial_star = best_trial_by_abs_dt(cache);
end

if isempty(trial_star) || ~trial_star.valid
    error('main_sweep_theta2:NoValidTrial', ...
        'Failed to obtain a valid trial for selecting the kinking angle.');
end

theta_star_deg = trial_star.theta_deg;
theta_star     = theta_star_deg*pi/180;

sig_star  = trial_star.sig_cr;
dt_star   = trial_star.dt_mouth;
ell_star  = trial_star.ell_adv;
j_star    = trial_star.j_star;

Pmid_star = build_Pmid(Pphys, dLeg, theta_star);

if S.verbose
    if bracket.found
        if bracket.iL == bracket.iR
            fprintf('root accepted directly at theta = %8.4f deg\n', bracket.thetaL_deg);
        else
            fprintf('bracket found after %d expansion(s): [%8.4f, %8.4f] deg\n', ...
                bracket.nExpand, bracket.thetaL_deg, bracket.thetaR_deg);
            fprintf('final bracket size  = %.6f deg\n', ...
                abs(bracket.thetaR_deg - bracket.thetaL_deg));
        end
        fprintf('refinement iterations = %d\n', trial_star.nRefine);
    else
        fprintf('no sign-changing bracket found; using minimum |dt_mouth| fallback\n');
    end

    dt_norm_star = nan;
    if isfield(C0,'Dmax') && ~isempty(C0.Dmax) && numel(C0.Dmax) >= 1 && abs(C0.Dmax(1)) > 0
        dt_norm_star = dt_star / C0.Dmax(1);
    end

    fprintf(['theta*            = %8.4f deg | beta = %.6g | sig = %.6g' ...
             ' | dt_mouth = % .3e | dt/dnmax = % .3e\n'], ...
        theta_star_deg, ell_star, sig_star, dt_star, dt_norm_star);

    fprintf('objective         = |dt_mouth| = %.3e\n', abs(dt_star));
end

% ---------------- optional exact re-solve at theta_star ----------------
out = struct();
if S.wantStarSol
    x0_star = trial_star.x;
    out.sol_star = solve_one_theta(Pphys, theta_star, C0, dLeg, x0_star, solverDisplay);
end

% ---------------- pack output ----------------
out.theta_star_deg = theta_star_deg;
out.theta_star     = theta_star;
out.sig_star       = sig_star;
out.dt_star        = dt_star;
out.ell_adv_star   = ell_star;
out.j_star         = j_star;
out.Pmid_star      = Pmid_star;

out.history   = cache_to_history_struct(cache);
out.bracket   = bracket;
out.predictor = pred;
out.cfg_used  = C0;

end


% ======================================================================
function pred = make_theta_predictor(theta_prev_deg, theta_prevprev_deg, theta_default_deg, beta, theta_min_deg, theta_max_deg)

if isempty(theta_prev_deg) || ~isfinite(theta_prev_deg)
    theta_pred_deg = theta_default_deg;
    trend_deg      = 0;
elseif isempty(theta_prevprev_deg) || ~isfinite(theta_prevprev_deg)
    theta_pred_deg = theta_prev_deg;
    trend_deg      = 0;
else
    trend_deg      = theta_prev_deg - theta_prevprev_deg;
    theta_pred_deg = theta_prev_deg + beta*trend_deg;
end

theta_pred_deg = min(max(theta_pred_deg, theta_min_deg), theta_max_deg);

pred = struct();
pred.theta_pred_deg = theta_pred_deg;
pred.trend_deg      = trend_deg;
pred.beta           = beta;
end


% ======================================================================
function [cache, idx] = eval_trial_cached(cache, theta_deg, Pphys, C0, dLeg, x0, solverDisplay, useWarmStart, verbose)

idx = find(abs(cache.theta_deg - theta_deg) < 1e-12, 1, 'first');
if ~isempty(idx)
    return;
end

if nargin < 6
    x0 = [];
end

% nearest successful warm start from cache, unless x0 was explicitly supplied
if isempty(x0) && useWarmStart
    x0 = nearest_successful_x(cache, theta_deg);
end

theta = theta_deg*pi/180;
sol = solve_one_theta(Pphys, theta, C0, dLeg, x0, solverDisplay);

trial = make_trial_struct(theta_deg, sol, dLeg);

cache.theta_deg(end+1,1) = theta_deg;
cache.trials{end+1,1}    = trial;

[cache.theta_deg, I] = sort(cache.theta_deg);
cache.trials = cache.trials(I);

idx = find(abs(cache.theta_deg - theta_deg) < 1e-12, 1, 'first');

if verbose
    if trial.valid
        if isfield(trial,'sol') && isstruct(trial.sol) && ...
                isfield(trial.sol,'cz') && isstruct(trial.sol.cz) && ...
                isfield(trial.sol.cz,'Dmax') && numel(trial.sol.cz.Dmax) >= 1 && ...
                abs(trial.sol.cz.Dmax(1)) > 0

        end

        fprintf(['  eval theta = %8.4f deg | flag = %d | it = %2d | sig = %.6g' ...
         ' | dt = % .3e | beta/delta = %.3f\n'], ...
        theta_deg, trial.flag, trial.iters, ...
        trial.sig_cr, trial.dt_mouth, trial.ell_adv / dLeg);
    else
        fprintf('  eval theta = %8.4f deg | flag = %d | FAILED\n', ...
            theta_deg, trial.flag);
    end
end
end


% ======================================================================
function x0 = nearest_successful_x(cache, theta_deg)

x0 = [];
if isempty(cache.trials)
    return;
end

bestDist = inf;
for k = 1:numel(cache.trials)
    tk = cache.trials{k};
    if tk.valid && ~isempty(tk.x)
        d = abs(tk.theta_deg - theta_deg);
        if d < bestDist
            bestDist = d;
            x0 = tk.x;
        end
    end
end
end


% ======================================================================
function trial = make_trial_struct(theta_deg, sol, dLeg)

trial = struct();
trial.theta_deg = theta_deg;
trial.theta     = theta_deg*pi/180;

trial.flag   = sol.flag;
trial.iters  = sol.iters;
trial.Fnorm  = sol.Fnorm;

trial.valid  = (sol.flag > 0);
trial.x      = [];
trial.sol    = sol;

trial.sig_cr   = nan;
trial.dt_mouth = nan;
trial.ell_adv  = nan;
trial.j_star   = nan;

trial.nRefine = 0;

if ~trial.valid
    return;
end

u     = sol.u;
cohes = sol.cohes;
cz    = sol.cz;

trial.sig_cr = sol.sig_cr;

% mouth opening (station 1 = mouth)
[~, Dt_m] = coh_opening_at_station(u, cohes, cz.Q, 1);
trial.dt_mouth = Dt_m;

% xi_n along stations, mouth -> tip
xi_n = compute_xi_n(u, cohes, cz.Q, cz);

% activated prefix
trial.j_star = active_zone_end_index(xi_n, cz.a1);
trial.ell_adv = crack_increment_from_index(trial.j_star, cohes, dLeg);

trial.x = sol.x;
end


% ======================================================================
function [bracket, cache] = find_bracket_adaptive(cache, Pphys, C0, dLeg, S, solverDisplay)

bracket = struct();
bracket.found      = false;
bracket.thetaL_deg = nan;
bracket.thetaR_deg = nan;
bracket.iL         = nan;
bracket.iR         = nan;
bracket.nExpand    = 0;

% first see whether the initial evaluations already contain a sign change
pair = find_first_sign_change_pair(cache);
if ~isempty(pair)
    bracket = finalize_bracket(pair, bracket, cache);
    return;
end

% decide preferred expansion direction from predictor trend or edge |dt|
dir = preferred_expansion_direction(cache, S);

step = S.expand_step_deg;

for k = 1:S.max_expand_iter
    bracket.nExpand = k;

    if dir < 0
        theta_try_deg = min(cache.theta_deg) - step;
        theta_try_deg = max(theta_try_deg, S.theta_min_deg);
    else
        theta_try_deg = max(cache.theta_deg) + step;
        theta_try_deg = min(theta_try_deg, S.theta_max_deg);
    end

    % if bound reached and point already exists, try opposite direction
    if any(abs(cache.theta_deg - theta_try_deg) < 1e-12)
        dir = -dir;
        if dir < 0
            theta_try_deg = min(cache.theta_deg) - step;
            theta_try_deg = max(theta_try_deg, S.theta_min_deg);
        else
            theta_try_deg = max(cache.theta_deg) + step;
            theta_try_deg = min(theta_try_deg, S.theta_max_deg);
        end
    end

    if any(abs(cache.theta_deg - theta_try_deg) < 1e-12)
        break;
    end

    [cache, ~] = eval_trial_cached(cache, theta_try_deg, Pphys, C0, dLeg, [], solverDisplay, S.useWarmStart, S.verbose);

    pair = find_first_sign_change_pair(cache);
    if ~isempty(pair)
        bracket = finalize_bracket(pair, bracket, cache);
        return;
    end

    step = step * S.expand_growth;
end
end


% ======================================================================
function dir = preferred_expansion_direction(cache, S)
% dir = -1 => expand to smaller angles first
% dir = +1 => expand to larger angles first

% 1) use history trend if available
if isfield(S,'theta_prev_deg') && isfield(S,'theta_prevprev_deg') && ...
        ~isempty(S.theta_prev_deg) && ~isempty(S.theta_prevprev_deg) && ...
        isfinite(S.theta_prev_deg) && isfinite(S.theta_prevprev_deg)

    tr = S.theta_prev_deg - S.theta_prevprev_deg;
    if tr < 0
        dir = -1;
        return;
    elseif tr > 0
        dir = +1;
        return;
    end
end

% 2) otherwise compare edge |dt|
valid = false(numel(cache.trials),1);
dt    = nan(numel(cache.trials),1);
for k = 1:numel(cache.trials)
    valid(k) = cache.trials{k}.valid;
    dt(k)    = cache.trials{k}.dt_mouth;
end

idxValid = find(valid);
if numel(idxValid) >= 2
    iL = idxValid(1);
    iR = idxValid(end);
    if abs(dt(iL)) <= abs(dt(iR))
        dir = -1;
    else
        dir = +1;
    end
else
    dir = -1;
end
end


% ======================================================================
function pair = find_first_sign_change_pair(cache)

pair = [];

n = numel(cache.trials);
if n < 2
    return;
end

for i = 1:(n-1)
    t1 = cache.trials{i};
    t2 = cache.trials{i+1};

    if ~(t1.valid && t2.valid)
        continue;
    end

    y1 = t1.dt_mouth;
    y2 = t2.dt_mouth;

    if ~isfinite(y1) || ~isfinite(y2)
        continue;
    end

    if y1 == 0
        pair = [i i];
        return;
    end
    if y2 == 0
        pair = [i+1 i+1];
        return;
    end

    if sign(y1) * sign(y2) < 0
        pair = [i i+1];
        return;
    end
end
end


% ======================================================================
function bracket = finalize_bracket(pair, bracket, cache)

bracket.found      = true;
bracket.iL         = pair(1);
bracket.iR         = pair(2);
bracket.thetaL_deg = cache.trials{pair(1)}.theta_deg;
bracket.thetaR_deg = cache.trials{pair(2)}.theta_deg;

end


% ======================================================================
function [trial_star, cache] = refine_bracketed_root(bracket, cache, Pphys, C0, dLeg, S, solverDisplay)

% exact-zero at cached point
if bracket.iL == bracket.iR
    trial_star = cache.trials{bracket.iL};
    trial_star.nRefine = 0;
    return;
end

iL = bracket.iL;
iR = bracket.iR;

trialL = cache.trials{iL};
trialR = cache.trials{iR};

thetaL_deg = trialL.theta_deg;
thetaR_deg = trialR.theta_deg;
yL = trialL.dt_mouth;
yR = trialR.dt_mouth;

nRef = 0;

while nRef < S.max_refine_iter
    nRef = nRef + 1;

    % stop if bracket already small
    if abs(thetaR_deg - thetaL_deg) <= S.theta_tol_deg
        break;
    end

    % regula-falsi step, safeguarded
    if isfinite(yL) && isfinite(yR) && abs(yR - yL) > eps
        thetaM_deg = thetaL_deg - yL*(thetaR_deg - thetaL_deg)/(yR - yL);
    else
        thetaM_deg = 0.5*(thetaL_deg + thetaR_deg);
    end

    % safeguard to interval interior
    thetaM_deg = max(thetaL_deg + 0.1*(thetaR_deg-thetaL_deg), ...
                 min(thetaR_deg - 0.1*(thetaR_deg-thetaL_deg), thetaM_deg));

    [cache, iM] = eval_trial_cached(cache, thetaM_deg, Pphys, C0, dLeg, [], solverDisplay, S.useWarmStart, S.verbose);
    trialM = cache.trials{iM};

    if ~trialM.valid
        % fallback to midpoint if regula-falsi landed on trouble
        thetaM_deg = 0.5*(thetaL_deg + thetaR_deg);
        [cache, iM] = eval_trial_cached(cache, thetaM_deg, Pphys, C0, dLeg, [], solverDisplay, S.useWarmStart, S.verbose);
        trialM = cache.trials{iM};

        if ~trialM.valid
            break;
        end
    end

    yM = trialM.dt_mouth;

    if abs(yM) <= S.dt_tol
        trial_star = trialM;
        trial_star.nRefine = nRef;
        return;
    end

    if sign(yL) * sign(yM) < 0
        thetaR_deg = trialM.theta_deg;
        yR = yM;
    else
        thetaL_deg = trialM.theta_deg;
        yL = yM;
    end
end

% choose best valid trial inside final bracket by min |dt|
trial_star = best_trial_in_interval_by_abs_dt(cache, thetaL_deg, thetaR_deg);
trial_star.nRefine = nRef;
end


% ======================================================================
function trial = best_trial_in_interval_by_abs_dt(cache, thetaL_deg, thetaR_deg)

trial = [];
bestVal = inf;

for k = 1:numel(cache.trials)
    tk = cache.trials{k};
    if ~tk.valid
        continue;
    end
    if tk.theta_deg < min(thetaL_deg, thetaR_deg) - 1e-12
        continue;
    end
    if tk.theta_deg > max(thetaL_deg, thetaR_deg) + 1e-12
        continue;
    end
    val = abs(tk.dt_mouth);
    if val < bestVal
        bestVal = val;
        trial = tk;
    end
end

if isempty(trial)
    trial = best_trial_by_abs_dt(cache);
end
end


% ======================================================================
function trial = best_trial_by_abs_dt(cache)

trial = [];
bestVal = inf;

for k = 1:numel(cache.trials)
    tk = cache.trials{k};
    if ~tk.valid || ~isfinite(tk.dt_mouth)
        continue;
    end
    val = abs(tk.dt_mouth);
    if val < bestVal
        bestVal = val;
        trial = tk;
    end
end
end


% ======================================================================
function H = cache_to_history_struct(cache)

n = numel(cache.trials);

H = struct();
H.theta_deg = nan(n,1);
H.flag      = nan(n,1);
H.sig_cr    = nan(n,1);
H.dt_mouth  = nan(n,1);
H.ell_adv   = nan(n,1);
H.beta      = H.ell_adv;
H.j_star    = nan(n,1);
H.iters     = nan(n,1);
H.Fnorm     = nan(n,1);

for k = 1:n
    tk = cache.trials{k};
    H.theta_deg(k) = tk.theta_deg;
    H.flag(k)      = tk.flag;
    H.sig_cr(k)    = tk.sig_cr;
    H.dt_mouth(k)  = tk.dt_mouth;
    H.ell_adv(k)   = tk.ell_adv;
    H.j_star(k)    = tk.j_star;
    H.iters(k)     = tk.iters;
    H.Fnorm(k)     = tk.Fnorm;
end
end


% ======================================================================
function Pmid = build_Pmid(Pphys, dLeg, theta)
tip  = Pphys(end,:);
Pmid = [Pphys; tip + dLeg*[cos(theta), sin(theta)]];
end


function Q = local_frame(Pmid)
v  = Pmid(end,:) - Pmid(end-1,:);
t  = v / max(norm(v), eps);
n  = [-t(2), t(1)];
Q = [t; n];
end


% ======================================================================
function xi_n = compute_xi_n(u, cohes, Q, cz)
mcoh  = size(cohes,1);
dnmax = cz.Dmax(1);
xi_n  = nan(mcoh,1);
for i = 1:mcoh
    [Dn, ~] = coh_opening_at_station(u, cohes, Q, i);
    xi_n(i) = Dn / dnmax;
end
end


function j_star = active_zone_end_index(xi_n, a1)
j_first_below = find(xi_n < a1, 1, 'first');
if isempty(j_first_below)
    j_star = numel(xi_n);
else
    j_star = max(1, j_first_below - 1);
end
end


function ell_adv = crack_increment_from_index(j_star, cohes, dLeg)
mcoh      = size(cohes,1);
ncoh_elem = (mcoh-1)/2;
ds = dLeg / (2*ncoh_elem);
ell_adv = (j_star - 1) * ds;
ell_adv = min(max(ell_adv, 0), dLeg);
end


% ======================================================================
function [Dn, Dt] = coh_opening_at_station(u, cohes, Q, i_station)
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


% ======================================================================
function sol = solve_one_theta(Pphys, theta, C0, dLeg, x0, solverDisplay)
%SOLVE_ONE_THETA Build geometry at a single theta and solve the critical state.

C = C0;

tip = Pphys(end,:);
C.Pmid = [Pphys; tip + dLeg*[cos(theta), sin(theta)]];
C.Q    = local_frame(C.Pmid);

% keep pencil width consistent with cohesive discretization
if ~isfield(C,'chw') || isempty(C.chw)
    C.chw = (dLeg / C.ncoh)/16;
end

[mesh, mat, quad, elod, Stif, cohes] = geom_pencil(C);

cz = struct('cohes',cohes, 'Q',C.Q, ...
    'Dmax',C.Dmax, 'sigmax',C.sigmax, ...
    'a1',C.a1, 'a2',C.a2, 'rphi',C.rphi);

ndof = 2*size(mesh.coord,1);

% Jacobian sparsity (bulk + CZ band + sigma)
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

x_scale = max(1e-9, (C.sig_guess/max(cz.sigmax(1),1))*C.B) * ones(ndof+1,1);

opts = optimoptions('fsolve', ...
    'Algorithm','trust-region-dogleg', ...
    'SpecifyObjectiveGradient', true, ...
    'JacobPattern', Jpat_aug, ...
    'ScaleProblem','jacobian', ...
    'TypicalX', x_scale, ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-12, ...
    'MaxIterations',100, ...
    'Display', char(solverDisplay));

if isempty(x0) || numel(x0) ~= (ndof+1)
    x0 = [zeros(ndof,1); C.sig_guess];
end

fun = @(xx) f1_crit3(xx, mesh, cz, Stif, elod, cohes);
[x, fval, flg, outfs] = fsolve(fun, x0, opts);

sol = struct();
sol.flag   = flg;
sol.iters  = outfs.iterations;
sol.Fnorm  = norm(fval);

sol.mesh   = mesh;
sol.mat    = mat;
sol.quad   = quad;
sol.elod   = elod;
sol.Stif   = Stif;
sol.cohes  = cohes;
sol.cz     = cz;

if flg > 0
    sol.u      = x(1:ndof);
    sol.sig_cr = x(end);
    sol.x      = x;
else
    sol.u      = nan(ndof,1);
    sol.sig_cr = nan;
    sol.x      = [];
end
end