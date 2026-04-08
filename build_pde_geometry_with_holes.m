function [mdl, dl, bt, gd, ns, sf] = build_pde_geometry_with_holes(BoundOuter, holeLoops)
%BUILD_PDE_GEOMETRY_WITH_HOLES Build PDE model from outer loop minus holes.
%
% Input
%   BoundOuter : N-by-2 outer boundary loop [x y]
%   holeLoops  : cell array of Nh-by-2 hole loops [x y]
%
% Output
%   mdl : PDE model created by createpde + geometryFromEdges
%   dl  : decomposed geometry matrix from decsg
%   bt  : Boolean table from decsg
%   gd  : geometry description matrix passed to decsg
%   ns  : names matrix passed to decsg
%   sf  : set formula passed to decsg
%
% Notes
%   - BoundOuter is standardized to counterclockwise orientation.
%   - Holes are standardized to clockwise orientation.
%   - If holeLoops is empty, this reproduces the single-region case.
%   - Each loop is treated as a polygon for decsg (type code 2).

    if nargin < 2 || isempty(holeLoops)
        holeLoops = {};
    end

    validate_loop(BoundOuter, 'BoundOuter');

    if ~iscell(holeLoops)
        error('build_pde_geometry_with_holes:BadInput', ...
            'holeLoops must be a cell array.');
    end

    % Clean and standardize outer loop
    BoundOuter = remove_duplicate_last_point(BoundOuter);
    BoundOuter = remove_consecutive_duplicates(BoundOuter);

    if size(BoundOuter,1) < 3
        error('build_pde_geometry_with_holes:DegenerateOuter', ...
            'BoundOuter must contain at least 3 unique vertices.');
    end

    % Standardize outer loop orientation: counterclockwise
    if polygon_signed_area(BoundOuter) < 0
        BoundOuter = flipud(BoundOuter);
    end

    % Clean and standardize hole loops: clockwise
    for k = 1:numel(holeLoops)
        validate_loop(holeLoops{k}, sprintf('holeLoops{%d}', k));

        H = remove_duplicate_last_point(holeLoops{k});
        H = remove_consecutive_duplicates(H);

        if size(H,1) < 3
            error('build_pde_geometry_with_holes:DegenerateHole', ...
                'holeLoops{%d} must contain at least 3 unique vertices.', k);
        end

        if polygon_signed_area(H) > 0
            H = flipud(H);
        end

        holeLoops{k} = H;
    end

    % Build geometry-description columns for decsg
    loops = [{BoundOuter}, holeLoops(:).'];
    nRegs = numel(loops);

    gdCols = cell(1, nRegs);
    names  = cell(1, nRegs);

    gdCols{1} = polygon_gd_column(BoundOuter);
    names{1}  = 'R1';

    for k = 1:numel(holeLoops)
        gdCols{k+1} = polygon_gd_column(holeLoops{k});
        names{k+1}  = sprintf('H%d', k);
    end

    gd = pack_gd_columns(gdCols);
    ns = char(names{:})';

    if isempty(holeLoops)
        sf = 'R1';
    else
        sf = ['R1' sprintf('-H%d', 1:numel(holeLoops))];
    end

    [dl, bt] = decsg(gd, sf, ns);

    mdl = createpde();
    geometryFromEdges(mdl, dl);
end


% -------------------------------------------------------------------------
function col = polygon_gd_column(P)
%POLYGON_GD_COLUMN Build one decsg polygon column for loop P.
%
% decsg polygon format:
%   [2;
%    Nv;
%    x1; x2; ...; xNv;
%    y1; y2; ...; yNv]
%
% Returned as a column vector.

    Nv = size(P,1);
    col = [2; Nv; P(:,1); P(:,2)];
end


% -------------------------------------------------------------------------
function gd = pack_gd_columns(cols)
%PACK_GD_COLUMNS Pad variable-length decsg columns and concatenate.

    lens = cellfun(@numel, cols);
    Lmax = max(lens);

    gd = zeros(Lmax, numel(cols));
    for j = 1:numel(cols)
        gd(1:lens(j), j) = cols{j};
    end
end


% -------------------------------------------------------------------------
function validate_loop(P, name)
%VALIDATE_LOOP Basic input validation for a polygon loop.

    if ~isnumeric(P) || size(P,2) ~= 2
        error('build_pde_geometry_with_holes:BadLoop', ...
            '%s must be an N-by-2 numeric array.', name);
    end
end


% -------------------------------------------------------------------------
function P = remove_duplicate_last_point(P)
    if size(P,1) >= 2
        if norm(P(end,:) - P(1,:), inf) == 0
            P(end,:) = [];
        end
    end
end


% -------------------------------------------------------------------------
function P = remove_consecutive_duplicates(P)
    keep = true(size(P,1),1);
    for i = 2:size(P,1)
        if norm(P(i,:) - P(i-1,:), inf) == 0
            keep(i) = false;
        end
    end
    P = P(keep,:);
end


% -------------------------------------------------------------------------
function A = polygon_signed_area(P)
    x = P(:,1);
    y = P(:,2);

    x2 = [x(2:end); x(1)];
    y2 = [y(2:end); y(1)];

    A = 0.5 * sum(x .* y2 - x2 .* y);
end