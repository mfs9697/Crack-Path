function LEFM = solve_next_angle_LEFM(Case)
%SOLVE_NEXT_ANGLE_LEFM  Compute LEFM kinking angle for one case.
%
% Inputs
%   Case   fully instantiated case struct from make_case_from_eta()
%
% Output
%   LEFM   result struct with LEFM kinking angle and optional diagnostics
%
% Method
%   - Build a one-leg straight-crack elastic problem using the current
%     non-perforated geometry and loading.
%   - Compute KI and KII from the circular J-integral with mode separation.
%   - Compute the LEFM kink angle using the MTS criterion.
%
% Notes
%   LEFM does not depend on the cohesive strengths directly. Therefore,
%   for fixed geometry and loading, the LEFM angle will normally be
%   constant across eta. We still compute it through this wrapper for
%   symmetry of the study workflow.

    LEFM = struct();

    LEFM.ok          = false;
    LEFM.msg         = '';
    LEFM.theta2_star = NaN;   % absolute angle of the predicted kinked leg
    LEFM.DeltaTheta  = NaN;   % relative kinking angle = theta2_star - theta1

    LEFM.theta_k_rel = NaN;   % signed local kink angle from MTS (rad)
    LEFM.theta_k     = NaN;   % alias of theta_k_rel for uniform downstream use
    LEFM.theta_k_deg = NaN;   % signed local kink angle (deg)

    LEFM.KI          = NaN;
    LEFM.KII         = NaN;

    LEFM.method      = 'MTS';
    LEFM.raw         = struct();

    try
        %% ================================================================
        %  build LEFM config in Crack Kinking format
        % ================================================================
        C = struct();

        % Geometry
        C.A    = Case.geometry.A;
        C.B    = Case.geometry.B;
        C.a    = Case.geometry.a;
        C.alf1 = Case.geometry.theta1;

        % Straight crack endpoints (one-leg problem)
        C.V0   = Case.geometry.P0;
        C.V1   = Case.geometry.P1;

        % Pencil channel width and mesh controls
        % Use the instantiated non-perforated mesh settings from Case.
        C.chw   = Case.mesh.chw;
        C.hmax  = Case.mesh.hmax;
        C.hgrad = Case.mesh.hgrad;

        % Edge-load smoothing
        % Keep this small and consistent with the old LEFM workflow.
        C.eps1 = 1e-9;

        % Material
        C.E2  = Case.material.E;
        C.nu  = Case.material.nu;
        C.G12 = Case.material.G;

        % Plane-strain constitutive matrix
        E  = C.E2;
        nu = C.nu;

        coef = E / ((1 + nu) * (1 - 2*nu));
        C.Dmat = coef * [ ...
            1-nu,   nu,          0; ...
            nu,     1-nu,        0; ...
            0,      0,   (1-2*nu)/2 ];

        % For optional visualization only
        C.plotMesh = false;

        % Cohesive-strength placeholders for plotting normalization only
        % (not physically needed by LEFM itself)
        C.sigmax = [Case.czm.sigmax_n, Case.czm.sigmax_t];

        %% ================================================================
        %  choose nominal reference load
        % ================================================================
        % The LEFM angle itself is independent of the nominal stress level,
        % but kinking_LEFM_1leg expects one. Use 1 MPa for clarity.
        sigma0 = 1.0;   % MPa

        %% ================================================================
        %  solve LEFM one-leg kinking problem
        % ================================================================
        [theta_k_rel, theta_k_deg, KI, KII] = kinking_LEFM_1leg(C, sigma0);

        theta2_abs = wrap_to_pi_local(C.alf1 - theta_k_rel);
        DeltaTheta = wrap_to_pi_local(theta2_abs - Case.theta1);

        %% ================================================================
        %  pack outputs
        % ================================================================
        LEFM.ok          = true;
        LEFM.msg         = 'success';

        LEFM.theta2_star = theta2_abs;
        LEFM.DeltaTheta  = DeltaTheta;

        LEFM.theta_k_rel = theta_k_rel;
        LEFM.theta_k     = theta_k_rel;   % alias for downstream consistency
        LEFM.theta_k_deg = theta_k_deg;

        LEFM.KI          = KI;
        LEFM.KII         = KII;

        LEFM.raw.C       = C;
        LEFM.raw.sigma0  = sigma0;

    catch ME
        LEFM.ok  = false;
        LEFM.msg = ME.message;
        LEFM.err = ME;
    end
end


% ========================================================================
function th = wrap_to_pi_local(th)
%WRAP_TO_PI_LOCAL  Wrap angle to (-pi, pi]
    th = mod(th + pi, 2*pi) - pi;
end