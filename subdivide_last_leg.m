function Pmid2 = subdivide_last_leg(Pmid, ncoh)
%SUBDIVIDE_LAST_LEG  Replace the last segment by ncoh segments.
% Pmid: (N x 2), N>=2
% Output: (N-1+ncoh x 2)

    N = size(Pmid,1);
    if N < 2
        error('Pmid must have at least 2 points.');
    end
    if ncoh < 1 || floor(ncoh) ~= ncoh
        error('ncoh must be a positive integer.');
    end

    P1 = Pmid(N-1,:);
    P2 = Pmid(N,:);

    t = linspace(0,1,ncoh+1).';
    seg = (1-t).*P1 + t.*P2;     % includes endpoints

    % Keep everything up to P1, then insert interior points, then P2
    Pmid2 = [Pmid(1:N-1,:); seg(2:end,:)];
end
