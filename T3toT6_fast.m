function [vert, tria6] = T3toT6_fast(vert, tria)
%T3TOT6_FAST  Upgrade T3 -> T6 using a global edge map.
% vert : [N x 2], tria : [Ne x 3]
% tria6: [Ne x 6] [v1 v2 v3 m12 m23 m31]

    ne = size(tria,1);
    n0 = size(vert,1);

    % all element edges (undirected)
    E12 = sort(tria(:,[1 2]),2);
    E23 = sort(tria(:,[2 3]),2);
    E31 = sort(tria(:,[3 1]),2);

    Eall = [E12; E23; E31];                % [3*Ne x 2]
    [Eu, ~, ic] = unique(Eall, 'rows');     % Eu: unique edges
    nu = size(Eu,1);

    % create mids for each unique edge
    mids = (n0+1 : n0+nu).';
    vert = [vert; 0.5*(vert(Eu(:,1),:) + vert(Eu(:,2),:))];

    % assign mids back to elements
    m12 = mids(ic(1:ne));
    m23 = mids(ic(ne+1:2*ne));
    m31 = mids(ic(2*ne+1:3*ne));

    tria6 = [tria, m12, m23, m31];
end
