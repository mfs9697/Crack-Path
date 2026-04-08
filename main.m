C = cfg_min_geom();

opts = struct('join', C.join, 'miter_limit', C.miter_limit);

Pmid = subdivide_last_leg(C.Pmid, C.ncoh);

[Bound, G] = build_domain_pencil_polyline(Pmid, C.A, C.B, C.chw, ...
    'join','miter', ...
    'miter_limit',6,...
    'tip', 'point',...    % <-- apex pencil
    'corner_tol',1e-10);            

figure(1); clf(1); hold on; axis equal; box on
plot([Bound(:,1); Bound(1,1)], [Bound(:,2); Bound(1,2)], 'k-');
plot(C.Pmid(:,1), C.Pmid(:,2), 'r.-', 'LineWidth', 1.5, 'MarkerSize', 14);
plot(G.up_chain(:,1), G.up_chain(:,2), 'b--');
plot(G.dn_chain(:,1), G.dn_chain(:,2), 'b--');
legend('Bound (domain boundary)','Midline','Offset +w','Offset -w');
title('Visible domain boundary with pencil cut');