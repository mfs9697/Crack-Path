function [f, M, Psi_hat] = tsl2(tD, cz)
%TSL2  Mixed-mode cohesive law (chi-based, conservative, symmetric core M).
%
% INPUT:
%   tD = [ xi_n , xi_t_signed ] where
%        xi_n = d_n/d_nmax (may be signed; xi_n<0 triggers contact penalty)
%        xi_t_signed = d_t/d_tmax (signed allowed; magnitude is used here)
%
% OUTPUT:
%   f      : [f_n; f_t] dimensionless traction magnitudes
%            f_n = sigma_n / sigma_n^max (>=0 on opening branch)
%            f_t = |sigma_t| / sigma_t^max  (>=0 ideally)
%            NOTE: shear sign is NOT included; apply sgn0(d_t) outside.
%
%   M      : 2x2 symmetric dimensionless "core" such that
%            K_loc = S * M * S,
%            S = diag( sqrt(sig_t^max/d_t^max), sqrt(sig_n^max/d_n^max) ).
%            (Ordering in M is [t, n] or [n, t]?  Here we use [n, t] to match f.)
%
%   Psi_hat: dimensionless mixed-mode potential (symmetric normalization),
%            Psi_hat = Psi / sqrt(chi_n*chi_t).
%
% THEORY (opening branch, shear magnitude):
%   Psi = chi_n F_n + chi_t F_t - r*sqrt(chi_n*chi_t)*F_n*F_t
%   sigma_n = sig_n^max * (1 - r*sqrt(chi_t/chi_n)*F_t) * t_n
%   |sigma_t|= sig_t^max * (1 - r*sqrt(chi_n/chi_t)*F_n) * t_t
%
% ADMISSIBILITY (nonnegative degradation for 0<=F<=c):
%   0 <= r <= rmax = min( (1/c_t)*sqrt(chi_n/chi_t), (1/c_n)*sqrt(chi_t/chi_n) ).
%
% NOTES:
% - g_base(xi,a1,a2,contactFlag) must return [F,t,g] for the given mode.
% - contactFlag=1 is used only for the normal mode (xi_n may be <0).
% - This function ignores the sign of d_t; caller applies sgn0(d_t).

% --- normalized openings
xi_n = tD(1);
xi_t = abs(tD(2));

% --- mode-wise shape parameters: a1 common, a2n/a2t may differ
a1 = cz.a1;

% Accept either cz.a2 = [a2n a2t] or explicit fields
if isfield(cz,'a2n') && isfield(cz,'a2t')
    a2n = cz.a2n;
    a2t = cz.a2t;
elseif isfield(cz,'a2') && numel(cz.a2) >= 2
    a2n = cz.a2(1);
    a2t = cz.a2(2);
else
    error('tsl2:MissingShape', ...
        'Need cz.a2=[a2n a2t] or fields cz.a2n and cz.a2t.');
end

% --- 1D base laws: (F,t,g). Normal includes contact penalty for xi_n<0.
[Fn, tn, gn] = g_base(xi_n, a1, a2n, 1);  % normal
[Ft, tt, gt] = g_base(xi_t, a1, a2t, 0);  % tangential (magnitude law)

% --- shape areas (dimensionless) for each mode
cn = (3 - 2*a1 + 3*a2n)/6;
ct = (3 - 2*a1 + 3*a2t)/6;

% --- determine energetic prefactors chi_n, chi_t
% Priority:
%   (1) if phi_n,phi_t provided: chi = phi/c
%   (2) else if sigmax & Dmax provided: chi = sigmax*Dmax
%   (3) legacy explicit fields: chi = sig_*max * D*max
if isfield(cz,'phi_n') && isfield(cz,'phi_t')
    phi_n = cz.phi_n;
    phi_t = cz.phi_t;
    chi_n = phi_n / cn;
    chi_t = phi_t / ct;

elseif isfield(cz,'sigmax') && isfield(cz,'Dmax')
    % cz.sigmax = [sig_nmax sig_tmax], cz.Dmax = [d_nmax d_tmax]
    chi_n = cz.sigmax(1) * cz.Dmax(1);
    chi_t = cz.sigmax(2) * cz.Dmax(2);

elseif all(isfield(cz,{'sig_nmax','sig_tmax','Dnmax','Dtmax'}))
    chi_n = cz.sig_nmax * cz.Dnmax;
    chi_t = cz.sig_tmax * cz.Dtmax;

else
    error('tsl2:MissingEnergyData', ...
        'Need (phi_n,phi_t) or (sigmax,Dmax) or legacy (sig_*max,D*max).');
end

chi_n = max(chi_n, realmin);
chi_t = max(chi_t, realmin);

% --- coupling parameter r
if isfield(cz,'r')
    r = cz.r;
elseif isfield(cz,'rphi')
    % allow legacy name, but interpret as r (NEW theory)
    r = cz.rphi;
else
    error('tsl2:MissingCoupling', 'Need coupling parameter cz.r (or legacy cz.rphi).');
end

% --- admissibility bound for r (nonnegative degradation on 0<=F<=c)
rmax_1 = (1/ct) * sqrt(chi_n/chi_t);
rmax_2 = (1/cn) * sqrt(chi_t/chi_n);
rmax   = min(rmax_1, rmax_2);

% Conservative check (you may change to warning if preferred)
tol = 1e-12;
if r < -tol
    error('tsl2:BadCoupling', 'Coupling r must be nonnegative. Got r=%g.', r);
end
if r > rmax + 1e-9
    error('tsl2:BadCoupling', ...
        'Coupling r=%g exceeds admissible bound rmax=%g (cn=%g, ct=%g).', r, rmax, cn, ct);
end

% --- dimensionless traction magnitudes
% Normalized by sigma_n^max and sigma_t^max, respectively.
Dn = 1 - r * sqrt(chi_t/chi_n) * Ft;  % normal degradation due to tangential separation
Dt = 1 - r * sqrt(chi_n/chi_t) * Fn;  % tangential degradation due to normal separation

fn = Dn * tn;   % sigma_n / sigma_n^max  (note: tn may be negative in compression branch)
ft = Dt * tt;   % |sigma_t| / sigma_t^max

f = [fn; ft];

% --- symmetric dimensionless core for K_loc = S*M*S
% IMPORTANT: M is returned in the same ordering as f = [n; t]
% so caller can map to full K_loc consistently.
%
% With ordering [n; t]:
%   M_nn = Dn * g_n
%   M_tt = Dt * g_t
%   M_nt = M_tn = - r * t_n * t_t
%
% If you prefer ordering [t; n] to match your paper, swap indices outside.
M = [ Dn * gn,      - r * tn * tt;
      - r * tn * tt, Dt * gt ];

% --- dimensionless potential (symmetric normalization)
% Psi_hat = Psi / sqrt(chi_n*chi_t)
Psi_dim = chi_n*Fn + chi_t*Ft - r*sqrt(chi_n*chi_t)*Fn*Ft;
Psi_hat = Psi_dim / sqrt(chi_n*chi_t);
end
