function Results = run_eta_kinking_study(C0)
%RUN_ETA_KINKING_STUDY  CZM vs LEFM kinking-angle study versus
% eta = E / sigmax_n for the non-perforated Crack Path branch.
%
% This version additionally computes an LEFM critical load using the
% quadratic energetic criterion via lefm_crit_load_vs_theta_quad().
%
% Usage:
%   Results = run_eta_kinking_study();
%   Results = run_eta_kinking_study(C0);

    %% ====================================================================
    %  baseline config
    % =====================================================================
    if nargin < 1 || isempty(C0)
        C0 = cfg_nonperforated();
    end

    eta_list = C0.study.eta_list(:);
    nEta     = numel(eta_list);

    if nEta == 0
        error('run_eta_kinking_study:EmptyEtaList', ...
            'C0.study.eta_list is empty.');
    end

    if isfield(C0, 'run') && isfield(C0.run, 'verbose')
        verbose = logical(C0.run.verbose);
    else
        verbose = true;
    end

    % LEFM quadratic critical-load settings
    if isfield(C0.study, 'lefm_sigma0') && ~isempty(C0.study.lefm_sigma0)
        lefm_sigma0 = C0.study.lefm_sigma0;
    else
        lefm_sigma0 = 1.0;  % MPa, nominal reference load
    end

    if isfield(C0.study, 'lefm_theta_grid_deg') && ~isempty(C0.study.lefm_theta_grid_deg)
        lefm_theta_grid_deg = C0.study.lefm_theta_grid_deg;
    else
        lefm_theta_grid_deg = (-90:1:90);
    end

    %% ====================================================================
    %  allocate result struct
    % =====================================================================
    Results = struct();

    Results.meta = struct();
    Results.meta.study_name   = C0.meta.name;
    Results.meta.description  = C0.meta.description;
    Results.meta.version      = C0.meta.version;
    Results.meta.timestamp    = char(datetime('now'));
    Results.meta.C0           = C0;

    Results.eta         = eta_list;
    Results.sigmax_n    = nan(nEta,1);
    Results.sigmax_t    = nan(nEta,1);

    Results.theta1      = nan(nEta,1);

    % CZM
    Results.theta2_CZM      = nan(nEta,1);
    Results.DeltaTheta_CZM  = nan(nEta,1);
    Results.sig_cr_CZM      = nan(nEta,1);
    Results.psi_hat_CZM     = nan(nEta,1);
    Results.dt_star_CZM     = nan(nEta,1);
    Results.ell_adv_CZM     = nan(nEta,1);
    Results.j_star_CZM      = nan(nEta,1);
    Results.beta_CZM        = nan(nEta,1);
    Results.beta_over_delta = nan(nEta,1);
    Results.ok_CZM          = false(nEta,1);
    Results.msg_CZM         = strings(nEta,1);
    Results.iters_CZM       = nan(nEta,1);

    % LEFM direction from existing solver
    Results.theta2_LEFM     = nan(nEta,1);
    Results.DeltaTheta_LEFM = nan(nEta,1);
    Results.KI_LEFM         = nan(nEta,1);
    Results.KII_LEFM        = nan(nEta,1);
    Results.ok_LEFM         = false(nEta,1);
    Results.msg_LEFM        = strings(nEta,1);

    % LEFM quadratic energetic critical load
    Results.theta_k_quad_LEFM_deg   = nan(nEta,1);   % local virtual kink angle, deg
    Results.theta2_quad_LEFM        = nan(nEta,1);   % absolute angle, rad
    Results.DeltaTheta_quad_LEFM    = nan(nEta,1);   % positive turn magnitude, rad
    Results.sig_cr_quad_LEFM        = nan(nEta,1);   % MPa
    Results.ok_quad_LEFM            = false(nEta,1);
    Results.msg_quad_LEFM           = strings(nEta,1);
    Results.theta_grid_quad_LEFM    = cell(nEta,1);
    Results.sigma_curve_quad_LEFM   = cell(nEta,1);

    % Differences
    Results.DeltaTheta_diff         = nan(nEta,1);   % CZM - LEFM(MTS)
    Results.DeltaTheta_diff_quad    = nan(nEta,1);   % CZM - LEFM(quadratic)
    Results.sig_cr_diff_CZM_quad    = nan(nEta,1);   % CZM - LEFM(quadratic)

    % Store raw case / solver outputs for traceability
    Results.Cases      = cell(nEta,1);
    Results.CZM        = cell(nEta,1);
    Results.LEFM       = cell(nEta,1);
    Results.LEFM_quad  = cell(nEta,1);

    %% ====================================================================
    %  LEFM caching logic
    % =====================================================================
    reuse_lefm = isfield(C0.study, 'compute_lefm_each') && ~C0.study.compute_lefm_each;
    lefm_cached = [];
    lefm_cached_ok = false;

    %% ====================================================================
    %  continuation state for CZM
    % =====================================================================
    Prev = [];

    %% ====================================================================
    %  main eta loop
    % =====================================================================
    if verbose
        fprintf('\n============================================================\n');
        fprintf('RUN_ETA_KINKING_STUDY\n');
        fprintf('Study: %s\n', C0.meta.name);
        fprintf('nEta = %d\n', nEta);
        fprintf('eta list = %s\n', mat2str(eta_list.'));
        fprintf('Crack angle (theta1) = %.1f deg\n', C0.geometry.theta1*180/pi);
        fprintf('============================================================\n');
    end

    for k = 1:nEta
        eta = eta_list(k);

        if verbose
            fprintf('\n------------------------------------------------------------\n');
            fprintf('ETA POINT %d / %d : eta = %.8g\n', k, nEta, eta);
            fprintf('------------------------------------------------------------\n');
        end

        %% ----------------------------------------------------------------
        %  build case
        % -----------------------------------------------------------------
        Case = make_case_from_eta(C0, eta);
        Case.id = k;

        Results.Cases{k}    = Case;
        Results.sigmax_n(k) = Case.czm.sigmax_n;
        Results.sigmax_t(k) = Case.czm.sigmax_t;
        Results.theta1(k)   = Case.theta1;

        %% ----------------------------------------------------------------
        %  LEFM solve (existing MTS-based direction)
        % -----------------------------------------------------------------
        if reuse_lefm && lefm_cached_ok
            LEFM = lefm_cached;
        else
            LEFM = solve_next_angle_LEFM(Case);

            if reuse_lefm && LEFM.ok
                lefm_cached    = LEFM;
                lefm_cached_ok = true;
            end
        end

        Results.LEFM{k}     = LEFM;
        Results.ok_LEFM(k)  = LEFM.ok;
        Results.msg_LEFM(k) = string(LEFM.msg);

        if LEFM.ok
            Results.theta2_LEFM(k)     = LEFM.theta2_star;
            Results.DeltaTheta_LEFM(k) = LEFM.DeltaTheta;

            if isfield(LEFM, 'KI')  && ~isempty(LEFM.KI),  Results.KI_LEFM(k)  = LEFM.KI;  end
            if isfield(LEFM, 'KII') && ~isempty(LEFM.KII), Results.KII_LEFM(k) = LEFM.KII; end

            %% ------------------------------------------------------------
            %  LEFM quadratic energetic critical load
            % -------------------------------------------------------------
            try
                KI0  = LEFM.KI;
                KII0 = LEFM.KII;

                [theta_grid_deg, sigma_quad_curve, theta_star_deg, sigma_star] = ...
                    lefm_crit_load_vs_theta_quad( ...
                    Case, KI0, KII0, ...
                    Case.czm.phi_n, Case.czm.phi_t, ...
                    lefm_sigma0, lefm_theta_grid_deg);

                % Use the quadratic criterion for the kink magnitude,
                % but orient it with the same turning side as LEFM(MTS).
                theta_star_rad_mag = deg2rad(abs(theta_star_deg));

                if isfield(LEFM, 'theta_k') && isfinite(LEFM.theta_k) && LEFM.theta_k ~= 0
                    kink_sign = sign(LEFM.theta_k);
                else
                    % fallback based on the standard small-angle LEFM sign
                    if isfinite(KII0) && KII0 ~= 0
                        kink_sign = -sign(KII0);
                    else
                        kink_sign = 1;
                    end
                end

                theta_k_quad = kink_sign * theta_star_rad_mag;
                theta2_quad  = Case.theta1 - theta_k_quad;
                Delta_quad   = abs(theta_k_quad);

                Results.theta_k_quad_LEFM_deg(k) = rad2deg(theta_k_quad);
                Results.theta2_quad_LEFM(k)      = theta2_quad;
                Results.DeltaTheta_quad_LEFM(k)  = Delta_quad;
                Results.sig_cr_quad_LEFM(k)      = sigma_star;
                Results.ok_quad_LEFM(k)          = isfinite(sigma_star) && isfinite(theta_star_deg);
                Results.msg_quad_LEFM(k)         = "ok";

                Results.theta_grid_quad_LEFM{k}  = theta_grid_deg;
                Results.sigma_curve_quad_LEFM{k} = sigma_quad_curve;

                Results.LEFM_quad{k} = struct( ...
                    'theta_deg',     theta_grid_deg, ...
                    'sigma_quad',    sigma_quad_curve, ...
                    'theta_star_deg',rad2deg(theta_k_quad), ...
                    'sigma_star',    sigma_star, ...
                    'theta2_star',   theta2_quad, ...
                    'DeltaTheta',    Delta_quad, ...
                    'sigma0',        lefm_sigma0);

            catch ME
                Results.ok_quad_LEFM(k)        = false;
                Results.msg_quad_LEFM(k)       = string(ME.message);
                Results.LEFM_quad{k}           = struct('error', ME.message);
            end

            %% ------------------------------------------------------------
            %  optional Dugdale-based delta update
            % -------------------------------------------------------------
            if isfield(C0.study, 'use_dugdale_delta') && C0.study.use_dugdale_delta && LEFM.ok

                KI  = LEFM.KI;
                KII = LEFM.KII;
                Keff = sqrt(KI^2 + KII^2);

                sigmax_n = Case.czm.sigmax_n;

                % user-specified fit constant and target beta/delta ratio
                C_D      = C0.study.dugdale_C;
                r_target = C0.study.beta_over_delta_target;

                rD = C_D * (Keff / sigmax_n)^2;
                delta_new = rD / r_target;

                % optional bounds for safety
                if isfield(C0.study, 'delta_min') && ~isempty(C0.study.delta_min)
                    delta_new = max(delta_new, C0.study.delta_min);
                end
                if isfield(C0.study, 'delta_max') && ~isempty(C0.study.delta_max)
                    delta_new = min(delta_new, C0.study.delta_max);
                end

                % update case
                Case.geometry.delta = delta_new;
                if isfield(Case.geometry, 'ds')
                    Case.geometry.ds = delta_new;
                end
                Case.delta = delta_new;

                Case.mesh.hmin = delta_new / Case.mesh.ncoh;
                Case.mesh.chw  = Case.mesh.hmin * C0.mesh.chw_factor;

                Case.derived.delta_over_a = delta_new / Case.geometry.a;

                Case.CP.delta = delta_new;
                Case.CP.hmin  = Case.mesh.hmin;
                Case.CP.chw   = Case.mesh.chw;

                % optional traceability
                Case.derived.Keff          = Keff;
                Case.derived.rD_est        = rD;
                Case.derived.delta_dugdale = delta_new;
            end
        else
            Results.ok_quad_LEFM(k)  = false;
            Results.msg_quad_LEFM(k) = "LEFM solve failed; quadratic post-processing skipped.";
        end

        %% ----------------------------------------------------------------
        %  CZM solve
        % -----------------------------------------------------------------
        CZM = solve_next_angle_CZM(Case, Prev);

        Results.CZM{k}     = CZM;
        Results.ok_CZM(k)  = CZM.ok;
        Results.msg_CZM(k) = string(CZM.msg);

        if CZM.ok
            Results.theta2_CZM(k)      = CZM.theta2_star;
            Results.DeltaTheta_CZM(k)  = CZM.DeltaTheta;
            Results.sig_cr_CZM(k)      = CZM.sig_cr;

            if isfield(CZM, 'psi_hat') && ~isempty(CZM.psi_hat)
                Results.psi_hat_CZM(k) = CZM.psi_hat;
            end
            if isfield(CZM, 'dt_star') && ~isempty(CZM.dt_star)
                Results.dt_star_CZM(k) = CZM.dt_star;
            end
            if isfield(CZM, 'ell_adv') && ~isempty(CZM.ell_adv)
                Results.ell_adv_CZM(k) = CZM.ell_adv;
            end
            if isfield(CZM, 'j_star') && ~isempty(CZM.j_star)
                Results.j_star_CZM(k) = CZM.j_star;
            end
            if isfield(CZM, 'beta') && ~isempty(CZM.beta)
                Results.beta_CZM(k) = CZM.beta;
            end
            if ~isfinite(Results.beta_CZM(k)) && isfinite(Results.ell_adv_CZM(k))
                Results.beta_CZM(k) = Results.ell_adv_CZM(k);
            end
            if isfield(CZM, 'beta_over_delta') && ~isempty(CZM.beta_over_delta)
                Results.beta_over_delta(k) = CZM.beta_over_delta;
            end
            if isfield(CZM, 'iters') && ~isempty(CZM.iters)
                Results.iters_CZM(k) = CZM.iters;
            end
        end

        %% ----------------------------------------------------------------
        %  differences
        % -----------------------------------------------------------------
        if CZM.ok && LEFM.ok
            Results.DeltaTheta_diff(k) = CZM.DeltaTheta - LEFM.DeltaTheta;
        end

        if CZM.ok && Results.ok_quad_LEFM(k)
            Results.DeltaTheta_diff_quad(k) = CZM.DeltaTheta - Results.DeltaTheta_quad_LEFM(k);
            Results.sig_cr_diff_CZM_quad(k) = CZM.sig_cr - Results.sig_cr_quad_LEFM(k);
        end

        %% ----------------------------------------------------------------
        %  update continuation predictor
        % -----------------------------------------------------------------
        if CZM.ok
            Prev = struct();
            Prev.ok          = true;
            Prev.theta2_star = CZM.theta2_star;

            if isfield(CZM, 'xsol') && ~isempty(CZM.xsol)
                Prev.xsol = CZM.xsol;
            else
                Prev.xsol = [];
            end

            if isfield(CZM, 'sig_cr') && ~isempty(CZM.sig_cr)
                Prev.sig_cr = CZM.sig_cr;
            else
                Prev.sig_cr = [];
            end
        else
            if ~isempty(Prev) && isfield(Prev, 'theta2_star') && isfinite(Prev.theta2_star)
                Prev.ok   = true;
                Prev.xsol = [];
            else
                Prev = [];
            end
        end

        %% ----------------------------------------------------------------
        %  logging
        % -----------------------------------------------------------------
        if verbose
            if LEFM.ok
                fprintf('LEFM(MTS)  : theta2 = %+10.6f deg | DeltaTheta = %+10.6f deg', ...
                    Results.theta2_LEFM(k)*180/pi, Results.DeltaTheta_LEFM(k)*180/pi);

                if isfinite(Results.KI_LEFM(k)) && isfinite(Results.KII_LEFM(k))
                    fprintf(' | KI = %.6g | KII = %.6g', Results.KI_LEFM(k), Results.KII_LEFM(k));
                end
                fprintf('\n');
            else
                fprintf('LEFM(MTS)  : FAILED | %s\n', LEFM.msg);
            end

            if Results.ok_quad_LEFM(k)
                theta2_deg   = rad2deg(Results.theta2_quad_LEFM(k));
                Delta_deg    = rad2deg(Results.DeltaTheta_quad_LEFM(k));

                fprintf('LEFM(quad) : theta_k = %+10.6f deg | theta2 = %+10.6f deg | DeltaTheta = %+10.6f deg | sig_cr = %.8g\n', ...
                    Results.theta_k_quad_LEFM_deg(k), ...
                    theta2_deg, ...
                    Delta_deg, ...
                    Results.sig_cr_quad_LEFM(k));
            else
                fprintf('LEFM(quad) : FAILED | %s\n', Results.msg_quad_LEFM(k));
            end

            if CZM.ok
                fprintf('CZM        : theta2 = %+10.6f deg | DeltaTheta = %+10.6f deg | sig_cr = %.8g', ...
                    Results.theta2_CZM(k)*180/pi, Results.DeltaTheta_CZM(k)*180/pi, Results.sig_cr_CZM(k));

                if isfinite(Results.beta_CZM(k))
                    fprintf(' | beta = %.6g', Results.beta_CZM(k));
                elseif isfinite(Results.ell_adv_CZM(k))
                    fprintf(' | beta = %.6g', Results.ell_adv_CZM(k));
                end

                if isfinite(Results.beta_over_delta(k))
                    fprintf(' | beta/delta = %.6g', Results.beta_over_delta(k));
                end

                if isfinite(Results.dt_star_CZM(k))
                    fprintf(' | dt_mouth = % .3e', Results.dt_star_CZM(k));
                end

                fprintf('\n');

                if isfinite(Results.DeltaTheta_diff(k))
                    fprintf('DIFF(MTS)  : DeltaTheta_CZM - DeltaTheta_LEFM(MTS)  = %+10.6f deg\n', ...
                        Results.DeltaTheta_diff(k)*180/pi);
                end

                if isfinite(Results.DeltaTheta_diff_quad(k))
                    fprintf('DIFF(quad) : DeltaTheta_CZM - DeltaTheta_LEFM(quad) = %+10.6f deg\n', ...
                        Results.DeltaTheta_diff_quad(k)*180/pi);
                end

                if isfinite(Results.sig_cr_diff_CZM_quad(k))
                    fprintf('LOAD DIFF  : sig_cr_CZM - sig_cr_LEFM(quad)         = %+10.6f MPa\n', ...
                        Results.sig_cr_diff_CZM_quad(k));
                end
            else
                fprintf('CZM        : FAILED | %s\n', CZM.msg);
            end
        end
    end

    %% ====================================================================
    %  summary
    % =====================================================================
    Results.summary = struct();
    Results.summary.nEta             = nEta;
    Results.summary.nCZM_success     = nnz(Results.ok_CZM);
    Results.summary.nLEFM_success    = nnz(Results.ok_LEFM);
    Results.summary.nLEFMquad_success= nnz(Results.ok_quad_LEFM);

    valid_diff = Results.DeltaTheta_diff(isfinite(Results.DeltaTheta_diff));
    if isempty(valid_diff)
        Results.summary.max_abs_diff_deg = NaN;
    else
        Results.summary.max_abs_diff_deg = max(abs(valid_diff))*180/pi;
    end

    valid_diff_quad = Results.DeltaTheta_diff_quad(isfinite(Results.DeltaTheta_diff_quad));
    if isempty(valid_diff_quad)
        Results.summary.max_abs_diff_quad_deg = NaN;
    else
        Results.summary.max_abs_diff_quad_deg = max(abs(valid_diff_quad))*180/pi;
    end

    valid_load_diff = Results.sig_cr_diff_CZM_quad(isfinite(Results.sig_cr_diff_CZM_quad));
    if isempty(valid_load_diff)
        Results.summary.max_abs_sig_cr_diff_CZM_quad = NaN;
    else
        Results.summary.max_abs_sig_cr_diff_CZM_quad = max(abs(valid_load_diff));
    end

    if verbose
        fprintf('\n============================================================\n');
        fprintf('Study complete.\n');
        fprintf('CZM success        : %d / %d\n', Results.summary.nCZM_success, nEta);
        fprintf('LEFM(MTS) success  : %d / %d\n', Results.summary.nLEFM_success, nEta);
        fprintf('LEFM(quad) success : %d / %d\n', Results.summary.nLEFMquad_success, nEta);
        fprintf('Max |DeltaTheta_CZM - DeltaTheta_LEFM(MTS)|  = %.6f deg\n', ...
            Results.summary.max_abs_diff_deg);
        fprintf('Max |DeltaTheta_CZM - DeltaTheta_LEFM(quad)| = %.6f deg\n', ...
            Results.summary.max_abs_diff_quad_deg);
        fprintf('Max |sig_cr_CZM - sig_cr_LEFM(quad)|         = %.6f MPa\n', ...
            Results.summary.max_abs_sig_cr_diff_CZM_quad);
        fprintf('============================================================\n');
    end

    %% ====================================================================
    %  save results if requested
    % =====================================================================
    if isfield(C0.study, 'save_results') && C0.study.save_results
        results_file = C0.study.results_file;
        save(results_file, 'Results');
        if verbose
            fprintf('Saved results to: %s\n', results_file);
        end
    end
end