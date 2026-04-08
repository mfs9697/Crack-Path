function h = fig_predictor_diagnostics(Path, varargin)
p = inputParser;
p.addParameter('kTurn', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('kLast', numel(Path.theta_star_deg), @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.parse(varargin{:});
S = p.Results;

k = (1:numel(Path.theta_star_deg)).';

h = figure('Color','w'); clf;
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

% --- theta* and predictor
nexttile; hold on; box on; grid on;
plot(k, Path.theta_star_deg, 'o-', 'LineWidth',1.4);
plot(k, Path.theta_pred_deg, 's--', 'LineWidth',1.2);
xline(S.kTurn, '--', 'LineWidth',1.2);
ylabel('deg');
legend({'$\theta^\star$','$\theta_{\mathrm{pred}}$'}, ...
    'Interpreter','latex', 'Location','best');
title('Predictor performance', 'Interpreter','latex');

% --- predictor error
nexttile; hold on; box on; grid on;
plot(k, Path.theta_pred_err_deg, 'o-', 'LineWidth',1.4);
yline(0, '-', 'LineWidth',1.0);
xline(S.kTurn, '--', 'LineWidth',1.2);
ylabel('$\theta^\star-\theta_{\mathrm{pred}}$', 'Interpreter','latex');

% --- total number of theta evaluations
nexttile; hold on; box on; grid on;
plot(k, Path.nThetaEval, 'o-', 'LineWidth',1.4);
xline(S.kTurn, '--', 'LineWidth',1.2);
xlabel('Propagation step $k$', 'Interpreter','latex');
ylabel('$n_{\theta,\mathrm{eval}}$', 'Interpreter','latex');
end