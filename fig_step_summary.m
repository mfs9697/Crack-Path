function h = fig_step_summary(step, Path, k, varargin)
%FIG_STEP_SUMMARY
% Three-panel figure:
%   (a) stress field
%   (b) critical load
%   (c) kinking angle

persistent hFig

p = inputParser;
p.addParameter('nstress', 2, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('stressFact', 10, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('figPosition', [100 100 1500 500], @(x)isnumeric(x)&&numel(x)==4);
p.addParameter('showTurnLine', false, @(b)islogical(b)&&isscalar(b));
p.parse(varargin{:});
S = p.Results;

if ~isfield(step,'sol_star') || ~isstruct(step.sol_star) || step.sol_star.flag <= 0
    error('fig_step_summary:NoStressSolution', ...
        'step.sol_star is missing or invalid for step %d.', k);
end

sol = step.sol_star;

% --- create or reuse figure ---
if isempty(hFig) || ~isvalid(hFig)
    hFig = figure('Color','w', 'Position', S.figPosition);
else
    figure(hFig);   % bring it to focus
end

clf(hFig);  % clear content but keep window

h = hFig;

TL = tiledlayout(h, 2, 2, 'Padding','compact', 'TileSpacing','compact');

% ============================================================
% (a) Stress field
% ============================================================
axA = nexttile(TL,[2 1]);
plot_stress_subfigure(axA, sol.mesh, sol.u, sol.cz, sol.mat, sol.quad, S.nstress, S.stressFact);
title(axA, sprintf('(a) Stress field at step $k=%d$', k), 'Interpreter','latex');

% 1:2 vertical-horizontal ratio
%axA.PlotBoxAspectRatio = [1 2 1];

% ============================================================
% (b) Critical load level
% ============================================================
axB = nexttile(TL);
hold(axB,'on'); box(axB,'on'); grid(axB,'on');

kk = (1:k).';
plot(axB, kk, Path.sig_star(1:k), 'o-', 'LineWidth', 1.5, 'MarkerSize', 5);
plot(axB, k, Path.sig_star(k), 'ro', 'MarkerFaceColor','r', 'MarkerSize', 7);

xlabel(axB, 'Propagation step $k$', 'Interpreter','latex');
ylabel(axB, '$\sigma_{\mathrm{cr}}\,[\mathrm{MPa}]$', 'Interpreter','latex');
title(axB, '(b) Critical load level', 'Interpreter','latex');

% ============================================================
% (c) Kinking angle
% ============================================================
axC = nexttile(TL,4);
hold(axC,'on'); box(axC,'on'); grid(axC,'on');

plot(axC, kk, Path.theta_star_deg(1:k), 'o-', 'LineWidth', 1.5, 'MarkerSize', 5);
plot(axC, k, Path.theta_star_deg(k), 'ro', 'MarkerFaceColor','r', 'MarkerSize', 7);

if S.showTurnLine && k >= 2
    [~, iMin] = min(Path.theta_star_deg(1:k));
    xline(axC, iMin, '--', 'LineWidth', 1.0);
end

xlabel(axC, 'Propagation step $k$', 'Interpreter','latex');
ylabel(axC, '$\theta_k^\star\,[\mathrm{deg}]$', 'Interpreter','latex');
title(axC, '(c) Kinking angle', 'Interpreter','latex');

end