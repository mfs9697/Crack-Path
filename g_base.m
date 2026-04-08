function [F, t, g] = g_base(xi, a1, a2, penalty)
% g_base
% -------------------------------------------------------------------------
% 1D cohesive base law returning a conservative triplet (F,t,g) such that
%   t = dF/dxi,   g = dt/dxi.
%
% Tension: trapezoidal-type law with saturation F(1) = c0,
%   c0 = \int_0^1 t(s) ds = (3 - 2*a1 + 3*a2)/6.
% Compression (if penalty=1): quadratic potential
%   F = (kappa/2) xi^2, t = kappa xi, g = kappa  for xi <= 0.
%
% INPUT:
%   xi      normalized separation (can be negative in compression)
%   a1,a2   TSL shape parameters (0<a1<a2<1)
%   penalty 1 enables contact penalty in compression, 0 disables
%
% OUTPUT:
%   F  dimensionless potential (NOT normalized to 1)
%   t  traction law (dimensionless)
%   g  tangent (dimensionless)

if nargin<4, penalty=0; end

% Simple linear penalty in compression: t = kappa*xi for xi<=0
if penalty
    kappa = (2/a1);
else
    kappa = 0;
end

% Shape parameter (= area under normalized traction curve on [0,1])
c0 = (3 - 2*a1 + 3*a2)/6;

% --- compression/contact branch: conservative quadratic potential
if xi <= 0
    F = 0.5*kappa*xi*xi;
    t = kappa*xi;
    g = kappa;
    return
end

a3 = (1 - a2)^3;

if xi < a1
    F = (xi^2/a1)*(1 - xi/(3*a1));
    t = (xi/a1)*(2 - xi/a1);
    g = (2/a1)*(1 - xi/a1);

elseif xi < a2
    F = xi - a1/3;
    t = 1;
    g = 0;

elseif xi <= 1
    xi2=xi*xi; xi3=xi2*xi; xi4=xi3*xi;
    A22=a2*a2; A23=A22*a2; A24=A23*a2;
    tm1=(xi4 - A24)/2 - (1 + a2)*(xi3 - A23) + 3*a2*(xi2 - A22) + (1 - 3*a2)*(xi - a2);
    F = a2 - a1/3 + tm1/a3;
    t = ((1 + 2*xi - 3*a2) * (1 - xi)^2) / a3;
    g = -6*(1 - xi)*(xi - a2) / a3;

else
    % fully broken: constant potential level, zero traction and tangent
    F = c0;
    t = 0;
    g = 0;
end
end
