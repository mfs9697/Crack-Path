function holeLoops = holes_to_loops(holes)
%HOLES_TO_LOOPS Convert hole specifications into polygonal boundary loops.
%
% Input
%   holes : empty, struct array, or cell array of structs.
%
% Supported hole specifications:
%   struct('type','circle',  'center',[xc yc], 'r',R, 'npoly',N)
%   struct('type','polygon', 'vertices',V)
%
% Output
%   holeLoops : cell array, each cell contains an Nh-by-2 loop [x y]
%
% Notes
%   - Each loop is returned as an open polygonal list of unique vertices:
%       [x1 y1;
%        x2 y2;
%        ...
%        xN yN]
%     without repeating the first point at the end.
%   - Orientation is standardized to clockwise for holes.
%   - Empty input returns {}.

    if nargin < 1 || isempty(holes)
        holeLoops = {};
        return;
    end

    % Normalize input to cell array of structs
    if isstruct(holes)
        holes = num2cell(holes);
    elseif ~iscell(holes)
        error('holes_to_loops:BadInput', ...
            'Input "holes" must be empty, a struct array, or a cell array of structs.');
    end

    holeLoops = cell(size(holes));

    for k = 1:numel(holes)
        hk = holes{k};

        if ~isstruct(hk) || ~isfield(hk, 'type')
            error('holes_to_loops:BadHoleSpec', ...
                'Each hole specification must be a struct with field "type".');
        end

        switch lower(strtrim(hk.type))
            case 'circle'
                validate_circle_spec(hk, k);

                c = hk.center(:).';
                R = hk.r;

                if isfield(hk, 'npoly') && ~isempty(hk.npoly)
                    npoly = hk.npoly;
                else
                    npoly = 80;  % default
                end

                H = circle_poly(c, R, npoly);

            case 'polygon'
                validate_polygon_spec(hk, k);
                H = hk.vertices;

            otherwise
                error('holes_to_loops:UnsupportedType', ...
                    'Unsupported hole type "%s" for hole %d.', hk.type, k);
        end

        H = remove_duplicate_last_point(H);
        H = remove_consecutive_duplicates(H);

        if size(H,1) < 3
            error('holes_to_loops:DegenerateLoop', ...
                'Hole %d has fewer than 3 unique vertices.', k);
        end

        % Standardize orientation for holes: clockwise
        if polygon_signed_area(H) > 0
            H = flipud(H);
        end

        holeLoops{k} = H;
    end
end


% -------------------------------------------------------------------------
function validate_circle_spec(h, k)
    req = {'center','r'};
    for i = 1:numel(req)
        if ~isfield(h, req{i})
            error('holes_to_loops:MissingField', ...
                'Circle hole %d is missing required field "%s".', k, req{i});
        end
    end

    if ~isnumeric(h.center) || numel(h.center) ~= 2
        error('holes_to_loops:BadCenter', ...
            'Circle hole %d: "center" must be a 1x2 or 2x1 numeric vector.', k);
    end

    if ~isscalar(h.r) || ~isnumeric(h.r) || h.r <= 0
        error('holes_to_loops:BadRadius', ...
            'Circle hole %d: "r" must be a positive numeric scalar.', k);
    end

    if isfield(h, 'npoly') && ~isempty(h.npoly)
        if ~isscalar(h.npoly) || h.npoly < 8 || round(h.npoly) ~= h.npoly
            error('holes_to_loops:BadNpoly', ...
                'Circle hole %d: "npoly" must be an integer >= 8.', k);
        end
    end
end


% -------------------------------------------------------------------------
function validate_polygon_spec(h, k)
    if ~isfield(h, 'vertices')
        error('holes_to_loops:MissingField', ...
            'Polygon hole %d is missing required field "vertices".', k);
    end

    V = h.vertices;
    if ~isnumeric(V) || size(V,2) ~= 2 || size(V,1) < 3
        error('holes_to_loops:BadVertices', ...
            'Polygon hole %d: "vertices" must be an N-by-2 numeric array, N >= 3.', k);
    end
end


% -------------------------------------------------------------------------
function H = circle_poly(c, R, n)
%CIRCLE_POLY Polygonal approximation of a circle.
% Returns an n-by-2 array of unique vertices, no repeated endpoint.

    t = linspace(0, 2*pi, n+1).';
    t(end) = [];

    H = [c(1) + R*cos(t), ...
         c(2) + R*sin(t)];
end


% -------------------------------------------------------------------------
function P = remove_duplicate_last_point(P)
%REMOVE_DUPLICATE_LAST_POINT Remove repeated closing point if present.

    if size(P,1) >= 2
        if norm(P(end,:) - P(1,:), inf) == 0
            P(end,:) = [];
        end
    end
end


% -------------------------------------------------------------------------
function P = remove_consecutive_duplicates(P)
%REMOVE_CONSECUTIVE_DUPLICATES Remove repeated neighboring vertices.

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
%POLYGON_SIGNED_AREA Signed area of polygon.
% Positive => CCW, negative => CW.

    x = P(:,1);
    y = P(:,2);

    x2 = [x(2:end); x(1)];
    y2 = [y(2:end); y(1)];

    A = 0.5 * sum(x .* y2 - x2 .* y);
end