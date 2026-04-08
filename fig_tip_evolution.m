function h = fig_tip_evolution(Path, varargin)
p = inputParser;
p.addParameter('kTurn', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('kLast', numel(Path.theta_star_deg), @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.parse(varargin{:});
S = p.Results;

k = (1:size(Path.tip_xy,1)).';

h = figure('Color','w'); clf;
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on; grid on;
plot(k, Path.tip_xy(:,1), 'o-', 'LineWidth',1.4);
xline(S.kTurn, '--', 'LineWidth',1.2);
ylabel('$x_{\mathrm{tip}}$', 'Interpreter','latex');
title('Evolution of crack-tip coordinates', 'Interpreter','latex');

nexttile; hold on; box on; grid on;
plot(k, Path.tip_xy(:,2), 'o-', 'LineWidth',1.4);
xline(S.kTurn, '--', 'LineWidth',1.2);
xlabel('Propagation step $k$', 'Interpreter','latex');
ylabel('$y_{\mathrm{tip}}$', 'Interpreter','latex');
end