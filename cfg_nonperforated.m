function C0 = cfg_nonperforated()
%CFG_NONPERFORATED  Baseline configuration for the non-perforated
% Crack Path study.
%
% This configuration is intended for:
%   (1) the IJF revision study:
%         dependence of the CZM kinking angle on eta = E/sigmax_n,
%         with comparison against the LEFM prediction;
%   (2) future TAFM work:
%         comparison of CZM and LEFM crack trajectories.
%
% Design principles:
%   * plain rectangular plate (no perforation);
%   * parameters aligned with the Crack Kinking paper-case as closely as
%     possible;
%   * trial cohesive-leg length delta is enlarged relative to the original
%     paper-case to avoid premature truncation of the active cohesive zone
%     in the eta-study.
%
% Units:
%   * stresses / moduli: MPa
%   * lengths: m
%   * fracture energies: J/m^2
%   * angles: rad

    C0 = struct();
    C0.eps1 = 1e-9;

    %% ====================================================================
    %  META
    % =====================================================================
    C0.meta = struct( ...
        'name',        'nonperforated_eta_study', ...
        'description', ['Non-perforated Crack Path branch aligned with ' ...
                        'Crack Kinking parameters; CZM vs LEFM kinking ' ...
                        'angle comparison versus eta = E/sigmax_n'], ...
        'version',     'v2_ck_aligned_large_delta');

    %% ====================================================================
    %  GEOMETRY
    % =====================================================================
    % Plate:
    %   x in [0, A], y in [-B, B]
    %
    % Crack:
    %   Leg 1 = physical crack, from P0 to P1
    %   Leg 2 = trial cohesive continuation, length delta, searched angle theta2
    %
    % Aligned with Crack Kinking:
    %   A = 0.30, B = 0.20, a = 0.05, theta1 = -20 deg
    %
    % Modified for eta-study:
    %   delta increased from 0.06*a to 0.12*a to better contain larger
    %   active cohesive zones at higher eta.

    C0.geometry = struct();

    C0.geometry.domain_type = 'rect_nonperforated';

    C0.geometry.A      = 0.10;
    C0.geometry.B      = 0.10;

    C0.geometry.P0     = [0.00, 0.00];
    C0.geometry.a      = 0.02;
    C0.geometry.theta1 = -25*pi/180;

    % Enlarged trial cohesive-leg length for applicability study
    C0.geometry.delta  = 0.08 * C0.geometry.a;

    % Reserved for future trajectory work
    C0.geometry.ds     = C0.geometry.delta;
    

    %% ====================================================================
    %  MATERIAL
    % =====================================================================
    % Crack Kinking paper-case:
    %   E = 4 GPa, nu = 0.3
    %
    % Here stored in MPa, consistent with the current MATLAB codebase.

    C0.material = struct();
    C0.eps1      = 1e-9;
    C0.sig_guess = 1.0;
    C0.plotMesh  = false;

    C0.material.E  = 4e3;   % MPa
    C0.material.nu = 0.30;

    %% ====================================================================
    %  COHESIVE MODEL
    % =====================================================================
    % Parameters aligned with the revised Crack Kinking manuscript:
    %   phi_n = 800 J/m^2
    %   phi_t = 1200 J/m^2
    %   a1    = 0.002
    %   a2n   = 0.90
    %   a2t   = 0.50
    %   rphi  = 0.40
    %
    % Reference strengths aligned with the paper-case:
    %   sigmax_n_ref = E/100
    %   sigmax_t_ref = G/100

    C0.czm = struct();

    % Fracture energies
    C0.czm.phi_n = 800;      % J/m^2
    C0.czm.phi_t = 1200;     % J/m^2

    % TSL shape parameters
    C0.czm.a1    = 0.002;
    C0.czm.a2n   = 0.90;
    C0.czm.a2t   = 0.50;

    % Mixed-mode energetic coupling
    C0.czm.rphi  = 0.40;

    % Descriptive only
    C0.czm.tsl_type = 'piecewise_smooth';

    % Reference strengths from Crack Kinking baseline
    Gref = C0.material.E / (2*(1 + C0.material.nu));

    C0.czm.sigmax_n_ref = C0.material.E / 100;   % 40 MPa
    C0.czm.sigmax_t_ref = Gref / 100;            % ~15.3846 MPa
    C0.czm.rho_sigma    = C0.czm.sigmax_t_ref / C0.czm.sigmax_n_ref;

    %% ====================================================================
    %  MESH
    % =====================================================================
    % Use refined paper-case settings as the default:
    %   ncoh = 80, hgrad = 1.15, B/hmax = 80
    %
    % Since delta is now larger, hmin = delta/ncoh is also larger unless
    % ncoh is increased. We keep ncoh = 80 initially and reassess after the
    % first eta sweep if needed.

    C0.mesh = struct();

    C0.mesh.ncoh       = 40;
    C0.mesh.hgrad      = 1.15;
    C0.mesh.hmax_ratio = 40;

    % Channel half-width rule:
    %   chw = hmin * chw_factor
    % Consistent with the current Crack Path convention.
    C0.mesh.chw_factor = 1/16;

    C0.mesh.order      = 'T6';

    %% ====================================================================
    %  ADAPTIVE CZM ANGLE SEARCH
    % =====================================================================
    % Search bounds in absolute angle theta2.
    %
    % Since theta1 = -20 deg and the selected absolute kink angle in the
    % Crack Kinking study lies near a few degrees, a moderate search window
    % around that regime is sufficient and more robust than a very wide one.

    C0.search = struct();

    C0.search.theta2_min = -15*pi/180;
    C0.search.theta2_max =  25*pi/180;

    % Initial local half-width around predictor
    C0.search.init_halfwidth = 1.0*pi/180;

    % Bracket expansion / refinement controls
    C0.search.expand_step_deg = 0.5;
    C0.search.expand_growth   = 2.0;
    C0.search.max_expand      = 8;

    C0.search.theta_tol_deg   = 0.01;
    C0.search.dt_tol          = 1e-12;
    C0.search.max_refine      = 20;

    % Use continuation in eta-space
    C0.search.use_predictor   = true;
    C0.search.use_warm_start  = true;

    % Default predictor if no previous eta-point exists
    C0.search.theta_default   = 2.0*pi/180;

    %% ====================================================================
    %  STUDY SETTINGS
    % =====================================================================
    % eta = E / sigmax_n
    %
    % Since sigmax_n_ref = E/100, the Crack Kinking paper-case corresponds
    % to eta = 100.
    %
    % We study both stronger and weaker cohesive interfaces relative to that
    % baseline.

    C0.study = struct();

    C0.study.eta_list = [50 75 100 125 150 200 250];

    C0.study.save_results      = true;
    C0.study.results_file      = 'eta_kinking_results.mat';

    % Optional diagnostics
    C0.study.compute_beta      = true;
    C0.study.compute_lefm_each = false;

    C0.study.use_dugdale_delta      = true;
    C0.study.dugdale_C              = 17.2*.8;
    C0.study.beta_over_delta_target = 0.9;

    % optional safety bounds
    C0.study.delta_min = 0.0006;
    C0.study.delta_max = 0.0020;

    %% ====================================================================
    %  LEFM SETTINGS
    % =====================================================================
    % LEFM angle will likely remain constant across eta for fixed geometry
    % and loading, but we keep the block explicit for future flexibility.

    C0.lefm = struct();

    C0.lefm.method       = 'MTS';
    C0.lefm.use_J_decomp = true;
    C0.lefm.theta_unit   = 'rad';

    %% ====================================================================
    %  INTEGRATION / RUN OPTIONS
    % =====================================================================
    C0.integration = struct();

    C0.integration.use_parent_path_functions = true;
    C0.integration.parent_root = '';

    C0.run = struct();
    C0.run.solverDisplay = 'off';
    C0.run.verbose       = true;

    %% ====================================================================
    %  DERIVED BASELINE QUANTITIES (REFERENCE ONLY)
    % =====================================================================
    % These are attached for logging / diagnostics and are not meant to
    % replace eta-dependent case construction.

    c_n = (3 - 2*C0.czm.a1 + 3*C0.czm.a2n) / 6;
    c_t = (3 - 2*C0.czm.a1 + 3*C0.czm.a2t) / 6;

    C0.derived = struct();

    C0.derived.G_ref   = Gref;
    C0.derived.c_n_ref = c_n;
    C0.derived.c_t_ref = c_t;

    % Convert MPa to Pa in the fracture-energy relation
    C0.derived.dmax_n_ref = C0.czm.phi_n / (c_n * C0.czm.sigmax_n_ref * 1e6);
    C0.derived.dmax_t_ref = C0.czm.phi_t / (c_t * C0.czm.sigmax_t_ref * 1e6);

    C0.derived.hmin_ref = C0.geometry.delta / C0.mesh.ncoh;
    C0.derived.hmax_ref = C0.geometry.B / C0.mesh.hmax_ratio;
    C0.derived.chw_ref  = C0.derived.hmin_ref * C0.mesh.chw_factor;

    % Baseline eta implied by the reference normal strength
    C0.derived.eta_ref = C0.material.E / C0.czm.sigmax_n_ref;

    % Convenience baseline points of the physical crack
    C0.derived.P1_ref = C0.geometry.P0 + ...
        C0.geometry.a * [cos(C0.geometry.theta1), sin(C0.geometry.theta1)];

end