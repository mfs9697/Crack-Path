function h = fig_trajectory_with_sig(Path, varargin)
%FIG_TRAJECTORY_WITH_SIG
% Two-panel figure:
%   (a) crack trajectory in the perforated domain
%   (b) critical load level sigma_cr vs propagation step
%
% Usage:
%   C0 = cfg_min_geom();
%   h = fig_trajectory_with_sig(Path, 'C0', C0, 'kTurn', 23);

p = inputParser;
p.addParameter('C0', [], @(x)isstruct(x) || isempty(x));
p.addParameter('kTurn', [], @(x)isempty(x) || (isnumeric(x)&&isscalar(x)&&x>=1));
p.addParameter('kLast', [], @(x)isempty(x) || (isnumeric(x)&&isscalar(x)&&x>=1));
p.addParameter('showHole', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('showTipMarkers', true, @(b)islogical(b)&&isscalar(b));
p.parse(varargin{:});
S = p.Results;

nK = numel(Path.theta_star_deg);
if isempty(S.kLast)
    S.kLast = nK;
end
if isempty(S.kTurn)
    [~, S.kTurn] = min(Path.theta_star_deg);
end

h = figure('Color','w'); clf;
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% ============================================================
% (a) Crack trajectory
% ============================================================
nexttile; hold on; axis equal; box on; grid on;

if ~isempty(S.C0)
    A = S.C0.A;
    B = S.C0.B;

    % outer boundary
    plot([0 A A 0 0], [-B -B B B -B], 'k-', 'LineWidth', 1.2);

    % hole(s)
    if S.showHole && isfield(S.C0,'holes') && ~isempty(S.C0.holes)
        holeLoops = holes_to_loops(S.C0.holes);
        for i = 1:numel(holeLoops)
            H = holeLoops{i};
            patch(H(:,1), H(:,2), [1 1 1], ...
                'EdgeColor', 'k', 'LineWidth', 1.2);
        end
    end
end

% initial crack
P0 = Path.Phist{1};
plot(P0(:,1), P0(:,2), 'k--', 'LineWidth', 1.2);

% full propagated crack path
Pfull = Path.Phist{end};
plot(Pfull(:,1), Pfull(:,2), 'r-', 'LineWidth', 2.0);

% selected markers
if S.showTipMarkers
    if S.kTurn <= size(Path.tip_xy,1) && all(isfinite(Path.tip_xy(S.kTurn,:)))
        plot(Path.tip_xy(S.kTurn,1), Path.tip_xy(S.kTurn,2), 'bo', ...
            'MarkerFaceColor','b', 'MarkerSize',7);
        text(Path.tip_xy(S.kTurn,1), Path.tip_xy(S.kTurn,2), ...
            sprintf('  k=%d', S.kTurn), 'Interpreter','none');
    end

    if S.kLast <= size(Path.tip_xy,1) && all(isfinite(Path.tip_xy(S.kLast,:)))
        plot(Path.tip_xy(S.kLast,1), Path.tip_xy(S.kLast,2), 'ms', ...
            'MarkerFaceColor','m', 'MarkerSize',7);
        text(Path.tip_xy(S.kLast,1), Path.tip_xy(S.kLast,2), ...
            sprintf('  k=%d', S.kLast), 'Interpreter','none');
    end
end

xlabel('$x$', 'Interpreter','latex');
ylabel('$y$', 'Interpreter','latex');
title('(a) Crack trajectory', 'Interpreter','latex');

% ============================================================
% (b) Critical load evolution
% ============================================================
nexttile; hold on; box on; grid on;

k = (1:numel(Path.sig_star)).';
plot(k, Path.sig_star, 'o-', 'LineWidth', 1.6, 'MarkerSize', 5);

xline(S.kTurn, '--', 'LineWidth', 1.2);
xline(S.kLast, ':', 'LineWidth', 1.2);

if S.kTurn <= numel(Path.sig_star) && isfinite(Path.sig_star(S.kTurn))
    plot(S.kTurn, Path.sig_star(S.kTurn), 'bo', 'MarkerFaceColor','b', 'MarkerSize',7);
end

if S.kLast <= numel(Path.sig_star) && isfinite(Path.sig_star(S.kLast))
    plot(S.kLast, Path.sig_star(S.kLast), 'ms', 'MarkerFaceColor','m', 'MarkerSize',7);
end

xlabel('Propagation step $k$', 'Interpreter','latex');
ylabel('$\sigma_{\mathrm{cr}}$', 'Interpreter','latex');
title('(b) Critical load level', 'Interpreter','latex');

end