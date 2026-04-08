function [theta_k, theta_deg, sigma_tt_val] = kink_angle_LEFM_MTS(KI, KII)
%KINK_ANGLE_LEFM_MTS  LEFM kinking angle for a 1-leg crack (MTS criterion).
%
%   [theta_k, theta_deg, sigma_tt_val] = kink_angle_LEFM_MTS(KI, KII)
%
%   theta_k     - rad, measured from the original crack plane
%   theta_deg   - deg
%   sigma_tt_val- circumferential stress (up to 1/sqrt(r) factor)
%                 at theta = theta_k

    % special cases
    if abs(KII) < 1e-12
        theta_k      = 0.0;
        theta_deg    = 0.0;
        sigma_tt_val = KI;
        return;
    end

    if abs(KI) < 1e-12
        % nearly pure mode II; use known ~70.5 deg as initial guess
        theta_guess = sign(KII) * deg2rad(70.5);
    else
        % small-angle estimate: tau_rt ~ KI*theta/2 + KII => theta ~ -2*KII/KI
        theta_guess = -2.0 * KII / KI;
        theta_guess = max(min(theta_guess, pi/2), -pi/2);
    end

    % stress components in polar coords (r factor omitted)
    function s_tt = sigma_tt(theta)
        c = cos(theta/2);
        s = sin(theta/2);
        s_tt = KI * c.^3 - 3.0 * KII .* s .* c.^2;
    end

    function t_rt = tau_rt(theta)
        c = cos(theta/2);
        s = sin(theta/2);
        t_rt = KI .* s .* c.^2 + KII .* c .* (1.0 - 3.0 * s.^2);
    end

    % solve tau_rt(theta) = 0 in a reasonable interval
    th_min = deg2rad(-80);
    th_max = deg2rad( 80);

    theta0 = min(max(theta_guess, th_min+1e-3), th_max-1e-3);
    tau_fun = @(th) tau_rt(th);

    dth = deg2rad(10);
    a = max(th_min, theta0 - dth);
    b = min(th_max, theta0 + dth);
    if sign(tau_fun(a)) * sign(tau_fun(b)) > 0
        a = th_min; b = th_max;
    end

    theta_k = fzero(tau_fun, [a, b]);
    sigma_tt_val = sigma_tt(theta_k);

    % if compressive, try opposite symmetric root (optional safeguard)
    if sigma_tt_val < 0
        theta_alt = -theta_k;
        if theta_alt > th_min && theta_alt < th_max
            sigma_tt_alt = sigma_tt(theta_alt);
            if sigma_tt_alt > sigma_tt_val
                theta_k      = theta_alt;
                sigma_tt_val = sigma_tt_alt;
            end
        end
    end

    theta_deg = rad2deg(theta_k);
end
