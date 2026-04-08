function h = fig_theta_evolution(Path, varargin)
p = inputParser;
p.addParameter('kTurn', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('kLast', numel(Path.theta_star_deg), @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.parse(varargin{:});
S = p.Results;

k = (1:numel(Path.theta_star_deg)).';

h = figure('Color','w'); clf; hold on; box on; grid on;
plot(k, Path.theta_star_deg, 'o-', 'LineWidth', 1.6, 'MarkerSize', 5);
xline(S.kTurn, '--', 'LineWidth', 1.2);
xline(S.kLast, ':', 'LineWidth', 1.2);

[thetaMin, iMin] = min(Path.theta_star_deg);
plot(iMin, thetaMin, 'ro', 'MarkerFaceColor','r');
text(iMin, thetaMin, sprintf('  min = %.2f^\\circ', thetaMin), 'Interpreter','tex');

xlabel('Propagation step $k$', 'Interpreter','latex');
ylabel('$\theta_k^\star\ [\mathrm{deg}]$', 'Interpreter','latex');
title('Evolution of the selected kinking angle', 'Interpreter','latex');
end