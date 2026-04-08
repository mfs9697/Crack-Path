function [theta_deg, sigma_quad, theta_star_deg, sigma_star] = ...
    lefm_crit_load_vs_theta_quad(C, KI0, KII0, GIc, GIIc, sigma0, theta_deg)
%LEFM_CRIT_LOAD_VS_THETA_QUAD
%   LEFM critical load versus virtual kink angle using the quadratic
%   energetic fracture criterion
%
%       GI/GIc + GII/GIIc = 1 .
%
%   This version:
%     1) controls the sign branch of theta automatically by default,
%     2) improves the precision of theta_star using coarse scan + fminbnd.
%
%   INPUTS:
%     C         : struct with fields
%                   .E or .E2  - Young's modulus (MPa)
%                   .nu        - Poisson ratio
%     KI0       : mode-I SIF of the original crack at reference load sigma0
%                 (MPa*sqrt(m))
%     KII0      : mode-II SIF of the original crack at reference load sigma0
%                 (MPa*sqrt(m))
%     GIc       : critical mode-I energy release rate (N/m)
%     GIIc      : critical mode-II energy release rate (N/m)
%     sigma0    : reference remote stress corresponding to KI0,KII0 (MPa)
%     theta_deg : optional coarse grid of virtual kink angles (deg),
%                 measured from the original crack plane.
%
%                 IMPORTANT:
%                 - If omitted, the function automatically chooses the
%                   physically relevant sign branch from KII0:
%                       KII0 > 0  --> search on negative branch [-90, 0]
%                       KII0 < 0  --> search on positive branch [0, 90]
%                 - If provided, the function respects the supplied grid.
%
%   OUTPUTS:
%     theta_deg     : column vector of tested coarse-grid angles (deg)
%     sigma_quad    : critical load for each coarse-grid angle (MPa)
%     theta_star_deg: refined angle at which sigma_quad is minimal (deg)
%     sigma_star    : minimal critical load (MPa)
%
%   NOTES:
%     - Local LEFM post-processing only.
%     - Uses Erdogan-Sih type small-kink SIF transformation.
%     - Energies are evaluated in plane strain:
%           Eeff = E/(1-nu^2).
%     - KI0 and KII0 must correspond to the same reference load sigma0.

    %% ---------------- Effective modulus ----------------
    if isfield(C, 'E') && ~isempty(C.E)
        E = C.E;
    elseif isfield(C, 'E2') && ~isempty(C.E2)
        E = C.E2;
    else
        error('lefm_crit_load_vs_theta_quad:MissingE', ...
            'C must contain either field E or E2.');
    end

    if ~isfield(C, 'nu') || isempty(C.nu)
        error('lefm_crit_load_vs_theta_quad:MissingNu', ...
            'C must contain field nu.');
    end

    Eeff = E / (1 - C.nu^2);   % MPa, plane strain

    %% ---------------- Default coarse grid with sign control ----------------
    if nargin < 7 || isempty(theta_deg)
        % Choose the physically relevant branch from the sign of KII0.
        % For the present sign convention and standard MTS small-angle logic:
        %   theta ~ -2*KII0/KI0
        % so the expected sign of theta is opposite to the sign of KII0
        if KII0 > 0
            theta_deg = (-90:1:0);
        elseif KII0 < 0
            theta_deg = (0:1:90);
        else
            % Pure mode I: symmetric case; include both sides
            theta_deg = (-90:1:90);
        end
    end

    theta_deg = theta_deg(:);   % ensure column

    if numel(theta_deg) < 3
        error('lefm_crit_load_vs_theta_quad:GridTooSmall', ...
            'theta_deg must contain at least 3 points.');
    end

    %% ---------------- Coarse evaluation ----------------
    sigma_quad = NaN(size(theta_deg));
    for i = 1:numel(theta_deg)
        sigma_quad(i) = sigma_quad_at_theta(theta_deg(i));
    end

    valid = isfinite(sigma_quad);
    if ~any(valid)
        theta_star_deg = NaN;
        sigma_star     = NaN;
        return;
    end

    [~, idx_valid_min] = min(sigma_quad(valid));
    idx_all_valid = find(valid);
    idx0 = idx_all_valid(idx_valid_min);

    %% ---------------- Local refinement with fminbnd ----------------
    % If the minimum is interior, refine in the neighboring bracket.
    % Otherwise, fall back to the coarse-grid minimum.
    if idx0 > 1 && idx0 < numel(theta_deg)
        a = theta_deg(idx0 - 1);
        b = theta_deg(idx0 + 1);

        obj = @(th_deg) sigma_quad_at_theta(th_deg);

        opts = optimset('TolX', 1e-6, 'Display', 'off');
        [theta_star_deg, sigma_star] = fminbnd(obj, a, b, opts);
    else
        theta_star_deg = theta_deg(idx0);
        sigma_star     = sigma_quad(idx0);
    end

    %% ---------------- Nested helper ----------------
    function sigma_val = sigma_quad_at_theta(th_deg)
        th = deg2rad(th_deg);

        c = cos(th/2);
        s = sin(th/2);

        % Erdogan-Sih small-kink transformation
        KIk  = KI0 * c^3 ...
             - 1.5 * KII0 * sin(th) * c;

        KIIk = KI0 * s * c^2 ...
             + KII0 * (cos(3*th/2) - 0.5 * c);

        % Energy release rates at sigma0
        % Units: MPa*m = 1e6 N/m
        GI0  = (KIk^2  / Eeff) * 1e6;   % N/m
        GII0 = (KIIk^2 / Eeff) * 1e6;   % N/m

        denom = GI0 / GIc + GII0 / GIIc;

        if denom > 0 && isfinite(denom)
            sigma_val = sigma0 / sqrt(denom);
        else
            sigma_val = NaN;
        end
    end
end