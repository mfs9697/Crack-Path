function [B, Det, dNdx] = BN_local(xi0, X)
% BN_LOCAL  T6 B-matrix at barycentric point (xi0) for element node coords X [6x2]
%   xi0 is 2x1 (standard triangle natural coords); xi3 = 1 - sum(xi0).
%
%   Inputs:
%     xi0 : [xi1; xi2]    (barycentric, xi3 = 1-xi1-xi2)
%     X   : [6 x 2]       (coordinates of T6 nodes)
%
%   Outputs:
%     B     : [3 x 12] strain-displacement matrix
%     Det   : determinant of Jacobian
%     dNdx  : [2 x 6] gradients of T6 shape functions w.r.t. x,y
%             (optional, for convenience)

    % barycentric coordinates
    xi  = [xi0(:); 1 - sum(xi0)];

    % Derivatives of shape functions wrt (xi1, xi2) for T6
    % L1 = xi(1), L2 = xi(2), L3 = xi(3)
    Nap = [ 4*xi(1)-1,      0,           1-4*xi(3),  4*xi(2),            -4*xi(2),            4*xi(3)-4*xi(1);
            0,              4*xi(2)-1,   1-4*xi(3),  4*xi(1),             4*xi(3)-4*xi(2),   -4*xi(1) ];

    % Jacobian of mapping: [dx/dxi1 dx/dxi2; dy/dxi1 dy/dxi2]
    dxdxi = Nap * X;            % [2 x 6] * [6 x 2] = [2 x 2]
    Det   = det(dxdxi);

    % Inverse Jacobian * Nap → gradients wrt x,y
    % N1 = [dN/dx; dN/dy]
    invJ = [ dxdxi(2,2), -dxdxi(1,2);
            -dxdxi(2,1),  dxdxi(1,1) ] / Det;

    dNdx = invJ * Nap;          % [2 x 2] * [2 x 6] = [2 x 6]

    % Build B-matrix
    nelvert = 6;
    eldf  = 2*nelvert; 
    inx = (2:2:eldf)'-1;        % [1,3,5,7,9,11]^T (ux DOFs)
    iny = inx+1;                % [2,4,6,8,10,12]^T (uy DOFs)

    B = zeros(3, eldf);
    B(1, inx) = dNdx(1,:);      % eps_xx = du/dx
    B(2, iny) = dNdx(2,:);      % eps_yy = dv/dy
    B(3, inx) = dNdx(2,:);      % gamma_xy = du/dy + dv/dx
    B(3, iny) = dNdx(1,:);

end

