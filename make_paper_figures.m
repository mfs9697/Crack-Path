function F = make_paper_figures(Path, varargin)
%MAKE_PAPER_FIGURES Build the main paper figures from Path.
%
% Usage:
%   C0 = cfg_min_geom();
%   F = make_paper_figures(Path, 'C0', C0, 'outDir', 'paper_figs');
%
% Optional name-value:
%   'C0'              : config from cfg_min_geom()
%   'outDir'          : directory for exported figures
%   'exportDPI'       : raster DPI for PNG
%   'kTurn'           : step index for maximum downward turning
%   'kLast'           : final step index
%   'showHole'        : true/false
%   'showTipMarkers'  : true/false

p = inputParser;
p.addParameter('C0', [], @(x)isstruct(x) || isempty(x));
p.addParameter('outDir', 'paper_figs', @(s)ischar(s)||isstring(s));
p.addParameter('exportDPI', 300, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('kTurn', [], @(x)isempty(x) || (isnumeric(x)&&isscalar(x)&&x>=1));
p.addParameter('kLast', [], @(x)isempty(x) || (isnumeric(x)&&isscalar(x)&&x>=1));
p.addParameter('showHole', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('showTipMarkers', true, @(b)islogical(b)&&isscalar(b));
p.parse(varargin{:});
S = p.Results;

if ~exist(S.outDir, 'dir')
    mkdir(S.outDir);
end

nK = numel(Path.theta_star_deg);
if isempty(S.kLast)
    S.kLast = nK;
end
if isempty(S.kTurn)
    [~, S.kTurn] = min(Path.theta_star_deg);
end

F = struct();

% 1) crack trajectory + critical load level
F.traj_sig = fig_trajectory_with_sig(Path, ...
    'C0', S.C0, ...
    'kTurn', S.kTurn, ...
    'kLast', S.kLast, ...
    'showHole', S.showHole, ...
    'showTipMarkers', S.showTipMarkers);
save_figure_pair(F.traj_sig, fullfile(S.outDir, 'fig01_trajectory_sig'), S.exportDPI);

% 2) kinking angle evolution
F.theta = fig_theta_evolution(Path, ...
    'kTurn', S.kTurn, ...
    'kLast', S.kLast);
save_figure_pair(F.theta, fullfile(S.outDir, 'fig02_theta_evolution'), S.exportDPI);

% 3) adaptive cohesive-length diagnostics
F.delta = fig_adaptive_delta(Path, ...
    'kTurn', S.kTurn, ...
    'kLast', S.kLast);
save_figure_pair(F.delta, fullfile(S.outDir, 'fig03_adaptive_delta'), S.exportDPI);

% 4) predictor / numerical diagnostics
F.pred = fig_predictor_diagnostics(Path, ...
    'kTurn', S.kTurn, ...
    'kLast', S.kLast);
save_figure_pair(F.pred, fullfile(S.outDir, 'fig04_predictor_diagnostics'), S.exportDPI);

fprintf('Saved paper figures to: %s\n', S.outDir);
fprintf('Suggested turning step: kTurn = %d\n', S.kTurn);
fprintf('Final step:             kLast = %d\n', S.kLast);
end


function save_figure_pair(figHandle, baseName, dpi)
savefig(figHandle, [baseName, '.fig']);
if exist('exportgraphics','file') == 2
    exportgraphics(figHandle, [baseName, '.png'], 'Resolution', dpi);
    exportgraphics(figHandle, [baseName, '.pdf'], 'ContentType', 'vector');
else
    print(figHandle, [baseName, '.png'], '-dpng', sprintf('-r%d', dpi));
end
end