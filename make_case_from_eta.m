function Case = make_case_from_eta(C0, eta)
%MAKE_CASE_FROM_ETA  Build one fully instantiated non-perforated case
% for a given modulus-to-strength ratio
%
%   eta = E / sigmax_n
%
% Inputs
%   C0   baseline configuration from cfg_nonperforated()
%   eta  scalar modulus-to-strength ratio
%
% Output
%   Case fully instantiated case struct for one eta-point
%
% Design rules
%   * both cohesive strengths vary proportionally
%   * fracture energies, TSL shape, coupling, geometry, and mesh policy
%     remain fixed
%   * all eta-dependent quantities are defined here, and nowhere else

    arguments
        C0   (1,1) struct
        eta  (1,1) double {mustBePositive}
    end

    %% ====================================================================
    %  initialize
    % =====================================================================
    Case = struct();
    Case.id  = [];
    Case.eta = eta;

    Case.eps1      = C0.eps1;
    Case.sig_guess = C0.sig_guess;
    Case.plotMesh  = C0.plotMesh;

    %% ====================================================================
    %  geometry
    % =====================================================================
    Case.geometry = C0.geometry;

    % Derived geometric points
    P0 = C0.geometry.P0;
    a  = C0.geometry.a;
    th1 = C0.geometry.theta1;

    P1 = P0 + a * [cos(th1), sin(th1)];
    P2_guess = P1 + C0.geometry.delta * [cos(th1), sin(th1)];

    Case.geometry.P1       = P1;
    Case.geometry.P2_guess = P2_guess;

    %% ====================================================================
    %  material
    % =====================================================================
    E  = C0.material.E;
    nu = C0.material.nu;
    G  = E / (2*(1 + nu));

    coef = E / ((1 + nu) * (1 - 2*nu));
    Dmat = coef * [ ...
    1-nu,   nu,          0; ...
    nu,     1-nu,        0; ...
    0,      0,   (1-2*nu)/2 ];

    Case.material = struct();
    Case.material.E  = E;
    Case.material.nu = nu;
    Case.material.G  = G;
    Case.material.Dmat = Dmat;

    %% ====================================================================
    %  cohesive model (eta-dependent strengths)
    % =====================================================================
    Case.czm = struct();

    % Fixed quantities copied from baseline
    Case.czm.phi_n  = C0.czm.phi_n;
    Case.czm.phi_t  = C0.czm.phi_t;

    Case.czm.a1     = C0.czm.a1;
    Case.czm.a2n    = C0.czm.a2n;
    Case.czm.a2t    = C0.czm.a2t;

    Case.czm.rphi   = C0.czm.rphi;
    Case.czm.rho_sigma = C0.czm.rho_sigma;
    Case.czm.tsl_type  = C0.czm.tsl_type;

    % Strengths derived from eta
    % eta = E / sigmax_n  => sigmax_n = E / eta
    sigmax_n = E / eta;
    sigmax_t = Case.czm.rho_sigma * sigmax_n;

    Case.czm.sigmax_n = sigmax_n;
    Case.czm.sigmax_t = sigmax_t;

    % Dimensionless areas under the normalized TSLs
    % c = (3 - 2*a1 + 3*a2)/6 for the adopted law
    c_n = (3 - 2*Case.czm.a1 + 3*Case.czm.a2n) / 6;
    c_t = (3 - 2*Case.czm.a1 + 3*Case.czm.a2t) / 6;

    Case.czm.c_n = c_n;
    Case.czm.c_t = c_t;

    % Critical separations from phi = sigmax * dmax * c
    % sigmax is in MPa, so convert to Pa in the denominator
    dmax_n = Case.czm.phi_n / (sigmax_n * c_n * 1e6);
    dmax_t = Case.czm.phi_t / (sigmax_t * c_t * 1e6);

    Case.czm.dmax_n = dmax_n;
    Case.czm.dmax_t = dmax_t;

    %% ====================================================================
    %  mesh
    % =====================================================================
    Case.mesh = C0.mesh;

    delta = C0.geometry.delta;
    B     = C0.geometry.B;
    ncoh  = C0.mesh.ncoh;

    hmin = delta / ncoh;
    hmax = B / C0.mesh.hmax_ratio;
    chw  = hmin * C0.mesh.chw_factor;

    Case.mesh.hmin = hmin;
    Case.mesh.hmax = hmax;
    Case.mesh.chw  = chw;

    %% ====================================================================
    %  search / LEFM blocks
    % =====================================================================
    Case.search = C0.search;
    Case.lefm   = C0.lefm;

    %% ====================================================================
    %  convenience shortcuts
    % =====================================================================
    Case.theta1 = Case.geometry.theta1;
    Case.delta  = Case.geometry.delta;

    Case.E  = Case.material.E;
    Case.nu = Case.material.nu;
    Case.G  = Case.material.G;

    %% ====================================================================
    %  seeds / predictor placeholders
    % =====================================================================
    Case.seed = struct();
    Case.seed.theta2_guess = [];
    Case.seed.xsol         = [];

    %% ====================================================================
    %  derived diagnostics / reference quantities
    % =====================================================================
    Case.derived = struct();

    % Relative trial-leg size
    Case.derived.delta_over_a = Case.geometry.delta / Case.geometry.a;

    % Baseline eta reference from cfg_nonperforated
    Case.derived.eta_ref = C0.derived.eta_ref;

    % Characteristic normal and tangential opening scales
    Case.derived.dmax = [dmax_n, dmax_t];

    % Strength vector in the order expected by most CZM routines
    Case.derived.sigmax = [sigmax_n, sigmax_t];

    % Fracture-energy vector
    Case.derived.phi = [Case.czm.phi_n, Case.czm.phi_t];

    % Shape-factor vector
    Case.derived.c = [c_n, c_t];

    %% ====================================================================
    %  Crack Path native configuration block
    % =====================================================================
    % This block is intended to minimize translation work in
    % solve_next_angle_CZM.m.
    %
    % IMPORTANT:
    %   The exact field names may need minor adjustment to match the current
    %   Crack Path main_sweep_theta2 / f1_crit3 conventions. But the key
    %   quantities are all assembled here.
    %
    % Convention:
    %   * stresses/moduli in MPa
    %   * lengths in m
    %   * Dmax in m
    %   * angles in rad

    CP = struct();
    CP.eps1      = Case.eps1;
    CP.sig_guess = Case.sig_guess;
    CP.plotMesh  = Case.plotMesh;
    CP.a2        = [Case.czm.a2n, Case.czm.a2t];
    

    % Geometry / loading
    CP.A     = Case.geometry.A;
    CP.B     = Case.geometry.B;
    CP.P0    = Case.geometry.P0;
    CP.a     = Case.geometry.a;
    CP.delta = Case.geometry.delta;

    % Leg-1 angle (keep both names for compatibility)
    CP.theta1 = Case.geometry.theta1;
    CP.alf1   = Case.geometry.theta1;

    % Material
    CP.E    = Case.material.E;
    CP.nu   = Case.material.nu;
    CP.G12  = Case.material.G;
    CP.E2   = Case.material.E;
    CP.Dmat = Case.material.Dmat;

    % Cohesive law
    CP.phi    = [Case.czm.phi_n, Case.czm.phi_t];
    CP.sigmax = [Case.czm.sigmax_n, Case.czm.sigmax_t];
    CP.sigmax_n = Case.czm.sigmax_n;
    CP.sigmax_t = Case.czm.sigmax_t;
    
    CP.Dmax   = [Case.czm.dmax_n,   Case.czm.dmax_t];
    CP.dmax_n = Case.czm.dmax_n;
    CP.dmax_t = Case.czm.dmax_t;
    
    CP.a1     = Case.czm.a1;
    CP.a2n    = Case.czm.a2n;
    CP.a2t    = Case.czm.a2t;
    CP.rphi   = Case.czm.rphi;

    % Shape-factor diagnostics (not necessarily used by solver directly)
    CP.c_n    = Case.czm.c_n;
    CP.c_t    = Case.czm.c_t;

    % Mesh
    CP.ncoh   = Case.mesh.ncoh;
    CP.hgrad  = Case.mesh.hgrad;
    CP.hmax   = Case.mesh.hmax;
    CP.hmin   = Case.mesh.hmin;
    CP.chw    = Case.mesh.chw;
    CP.order  = Case.mesh.order;

    % Search settings
    CP.theta2_min     = Case.search.theta2_min;
    CP.theta2_max     = Case.search.theta2_max;
    CP.init_halfwidth = Case.search.init_halfwidth;
    CP.expand_step_deg = Case.search.expand_step_deg;
    CP.expand_growth   = Case.search.expand_growth;
    CP.max_expand      = Case.search.max_expand;
    CP.theta_tol_deg   = Case.search.theta_tol_deg;
    CP.dt_tol          = Case.search.dt_tol;
    CP.max_refine      = Case.search.max_refine;
    CP.theta_default   = Case.search.theta_default;

    % Optional run flags
    CP.use_predictor   = Case.search.use_predictor;
    CP.use_warm_start  = Case.search.use_warm_start;

    % Eta study traceability
    CP.eta             = Case.eta;

    Case.CP = CP;

    %% ====================================================================
    %  validation checks
    % =====================================================================
    % Basic sanity checks that should fail early if the case is inconsistent.

    assert(Case.czm.sigmax_n > 0, 'sigmax_n must be positive.');
    assert(Case.czm.sigmax_t > 0, 'sigmax_t must be positive.');

    assert(Case.czm.dmax_n > 0, 'dmax_n must be positive.');
    assert(Case.czm.dmax_t > 0, 'dmax_t must be positive.');

    assert(Case.mesh.hmin > 0, 'hmin must be positive.');
    assert(Case.mesh.hmax > 0, 'hmax must be positive.');
    assert(Case.mesh.chw  > 0, 'chw must be positive.');

    assert(Case.search.theta2_min < Case.search.theta2_max, ...
        'theta2_min must be smaller than theta2_max.');

end