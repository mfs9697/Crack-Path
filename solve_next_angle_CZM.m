function CZM = solve_next_angle_CZM(Case, Prev)
%SOLVE_NEXT_ANGLE_CZM  Solve one CZM kinking-angle problem for a given eta.
%
% Wrapper around the Crack Path adaptive angle-search engine.
%
% Inputs
%   Case   fully instantiated case struct from make_case_from_eta()
%   Prev   previous eta-point result struct for predictor / warm start
%          (may be empty [])
%
% Output
%   CZM    result struct with selected angle, critical load, diagnostics,
%          and optional warm-start data for the next eta-point
%
% Notes
%   This wrapper assumes the adaptive Crack Path search engine can accept
%   the case-native config struct through a name-value pair:
%
%       main_sweep_theta2(Pphys, 'cfg', Case.CP, ...)
%
%   If the parent project does not yet support this, create a local sibling
%   function (e.g. main_sweep_theta2_case.m) with the same API and replace
%   the call below accordingly.

    if nargin < 2
        Prev = [];
    end

    %% ====================================================================
    %  initialize output
    % =====================================================================
    CZM = struct();

    CZM.ok         = false;
    CZM.msg        = '';
    CZM.theta2_star = NaN;
    CZM.DeltaTheta  = NaN;
    CZM.sig_cr      = NaN;
    CZM.dt_star     = NaN;
    CZM.ell_adv     = NaN;
    CZM.j_star      = NaN;
    CZM.psi_hat     = NaN;

    CZM.iters       = NaN;
    CZM.theta_bracket = [];
    CZM.history     = [];
    CZM.predictor   = [];
    CZM.Pmid_star   = [];
    CZM.sol_star    = [];
    CZM.xsol        = [];

    % Optional active cohesive-zone diagnostics
    CZM.beta        = NaN;
    CZM.beta_over_delta = NaN;

    %% ====================================================================
    %  build physical crack polyline Pphys (mouth -> current tip)
    % =====================================================================
    % For the non-perforated single-step study, the physical crack is just
    % Leg 1, from P0 to P1.
    %
    % Later, for trajectory work, this can be generalized to a multi-point
    % physical crack path.

    Pphys = [Case.geometry.P0; Case.geometry.P1];

    %% ====================================================================
    %  predictor / warm start from previous eta-point
    % =====================================================================
    theta_prev_deg     = [];
    theta_prevprev_deg = []; %#ok<NASGU> % reserved for future 2-step continuation
    x0_prev            = [];

    if ~isempty(Prev) && isstruct(Prev)
        if isfield(Prev, 'ok') && Prev.ok
            if isfield(Prev, 'theta2_star') && isfinite(Prev.theta2_star)
                theta_prev_deg = Prev.theta2_star * 180/pi;
            end
            if isfield(Prev, 'xsol') && ~isempty(Prev.xsol)
                x0_prev = Prev.xsol;
            end
        end
    end

    %% ====================================================================
    %  select search engine
    % =====================================================================
    % Preferred:
    %   call the parent Crack Path adaptive engine with Case.CP injected.
    %
    % If you instead create a local copied/adapted engine in the new folder,
    % replace the function handle below.

    %% ====================================================================
    %  call adaptive CZM search
    % =====================================================================
    try
        % ---- predictor window ----
        theta_default_deg = Case.search.theta_default * 180/pi;

        % If no previous point exists, center around default predictor.
        if isempty(theta_prev_deg)
            theta_prev_deg = theta_default_deg;
        end

        % Note:
        % main_sweep_theta2 currently supports asymmetric initial windows
        % through init_left_deg / init_right_deg.
        init_half_deg = Case.search.init_halfwidth * 180/pi;

        out = sweep_theta2( ...
                Pphys, ...
                'cfg',               Case.CP, ...
                'delta',             Case.delta, ...
                'theta_prev_deg',    theta_prev_deg, ...
                'theta_prevprev_deg', [], ...
                'theta_default_deg', theta_default_deg, ...
                'init_left_deg',     init_half_deg, ...
                'init_right_deg',    init_half_deg, ...
                'trendBias',         1.0, ...
                'theta_min_deg',     Case.search.theta2_min * 180/pi, ...
                'theta_max_deg',     Case.search.theta2_max * 180/pi, ...
                'expand_step_deg',   Case.search.expand_step_deg, ...
                'expand_growth',     Case.search.expand_growth, ...
                'max_expand_iter',   Case.search.max_expand, ...
                'theta_tol_deg',     Case.search.theta_tol_deg, ...
                'dt_tol',            Case.search.dt_tol, ...
                'max_refine_iter',   Case.search.max_refine, ...
                'solverDisplay',     'off', ...
                'useWarmStart',      Case.search.use_warm_start, ...
                'wantStarSol',       true, ...
                'verbose',           true);

        %% ================================================================
        %  pack outputs
        % ================================================================
        CZM.ok          = true;
        CZM.msg         = 'success';

        CZM.theta2_star = out.theta_star;
        CZM.DeltaTheta  = wrap_to_pi_local(out.theta_star - Case.theta1);

        CZM.sig_cr      = out.sig_star;

        CZM.dt_star     = out.dt_star;
        CZM.beta        = out.ell_adv_star;   % beta in the Crack Kinking paper
        CZM.ell_adv     = CZM.beta;           % backward-compatible alias
        CZM.beta_over_delta = CZM.beta / Case.delta;
        CZM.j_star      = out.j_star;

        CZM.Pmid_star   = out.Pmid_star;

        if isfield(out, 'sol_star')
            CZM.sol_star = out.sol_star;
            CZM.xsol     = out.sol_star.x;
        end

        if isfield(out, 'history')
            CZM.history = out.history;
        end
        if isfield(out, 'bracket')
            CZM.theta_bracket = out.bracket;
        end
        if isfield(out, 'predictor')
            CZM.predictor = out.predictor;
        end

        % Try to infer number of refinement steps if present
        if isfield(out, 'bracket') && isstruct(out.bracket)
            if isfield(out.bracket, 'nExpand')
                CZM.iters = out.bracket.nExpand;
            end
        end

        %% ================================================================
        %  optional postprocessing: active cohesive-zone length beta
        % ================================================================
        if isfield(Case, 'study') && isfield(Case.study, 'compute_beta') && Case.study.compute_beta
            % left here for future integration if Case carries the flag
        end

        % If available, extract beta from the converged solution
        if ~isempty(CZM.sol_star)
            try
                [beta, beta_over_delta] = extract_active_cz_length(CZM.sol_star, Case);
                CZM.beta = beta;
                CZM.beta_over_delta = beta_over_delta;
            catch
                % keep NaN if helper not ready yet
            end
        end

        %% ================================================================
        %  optional normalized mouth energy diagnostic
        % ================================================================
        % main_sweep_theta2 currently returns dt-based root information but
        % not psi_hat explicitly. Keep as NaN unless a local engine or
        % helper computes it.
        %
        % If later the local case-engine returns psi_hat, just map it here:
        %
        % if isfield(out,'psi_hat_star'), CZM.psi_hat = out.psi_hat_star; end

    catch ME
        CZM.ok   = false;
        CZM.msg  = ME.message;
        CZM.err  = ME;
    end
end


% ========================================================================
function th = wrap_to_pi_local(th)
% local wrap to (-pi, pi]
    th = mod(th + pi, 2*pi) - pi;
end