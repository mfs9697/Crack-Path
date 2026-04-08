function h = fig_adaptive_delta(Path, varargin)
p = inputParser;
p.addParameter('kTurn', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('kLast', numel(Path.theta_star_deg), @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.parse(varargin{:});
S = p.Results;

k = (1:numel(Path.ell_adv)).';

h = figure('Color','w'); clf;

tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on; grid on;
plot(k, Path.ell_adv, 'o-', 'LineWidth',1.5, 'MarkerSize',5);
plot(k, Path.delta_used, 's-', 'LineWidth',1.5, 'MarkerSize',5);
xline(S.kTurn, '--', 'LineWidth', 1.2);
legend({'$\ell_{\mathrm{adv}}$','$\delta$'}, 'Interpreter','latex', 'Location','best');
ylabel('Length', 'Interpreter','latex');
title('Adaptive cohesive-length diagnostics', 'Interpreter','latex');

nexttile; hold on; box on; grid on;
plot(k, Path.active_ratio, 'o-', 'LineWidth',1.5, 'MarkerSize',5);
yline(0.8, '--', 'LineWidth',1.0);
yline(0.9, '--', 'LineWidth',1.0);
xline(S.kTurn, '--', 'LineWidth', 1.2);
xlabel('Propagation step $k$', 'Interpreter','latex');
ylabel('$\ell_{\mathrm{adv}}/\delta$', 'Interpreter','latex');
end