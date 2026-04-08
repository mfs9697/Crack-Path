function [KI, KII] = SIF_LEFM_circle2(mesh, U, V, mat, rI)
%SIF_LEFM_CIRCLE  Compute KI, KII using circular J-integral with
%                 Ishikawa–Kitagawa–Okamura mode separation.
%
%   [KI, KII] = SIF_LEFM_circle(mesh, U, V, mat, rI)
%
%   Inputs:
%     mesh   : struct with fields
%              .coord    [nT6 x 2]  T6 node coordinates
%              .connect  [nElem x 6] T6 connectivity
%              .coord3   [nT3 x 2]  T3 parent coordinates (corner nodes)
%              .connect3 [nElem x 3] T3 parent connectivity
%     U      : displacement vector (2*nnod x 1), [ux1; uy1; ux2; ...]
%     V      : [nPts x 2] crack polyline, last point = crack tip
%     mat    : struct with .E, .nu, .Dmat  (plane strain stiffness)
%     rI     : contour radius (if empty -> 0.7 * total crack length)
%
%   Output:
%     KI, KII : stress intensity factors (plane strain)

    % ---------------------------------------------------------------------
    % 0.  GET MESH, MATERIAL
    % ---------------------------------------------------------------------
    coord    = mesh.coord;      % T6
    connect  = mesh.connect;    % T6
    coord3   = mesh.coord3;     % T3
    connect3 = mesh.connect3;   % T3

    E    = mat.E;
    nu   = mat.nu;
    Dmat = mat.Dmat;            % plane strain isotropic stiffness

    nelvert = 6;
    eldf    = 2*nelvert;
    iny     = (2:2:eldf)';      % [2,4,6,8,10,12]^T
    inx     = iny - 1;          % [1,3,5,7,9,11]^T

    % ---------------------------------------------------------------------
    % 1.  CRACK TIP, CRACK DIRECTION, DEFAULT RADIUS
    % ---------------------------------------------------------------------
    x_tip = V(end,:).';              % crack tip (2x1)
    a_tot = norm(V(end,:).'-V(1,:).');  % total crack length

    if nargin < 5 || isempty(rI)
        rI = 0.7 * a_tot;
    end

    % Crack tangent (global) from last segment
    p_tip  = V(end,:).';
    p_prev = V(end-1,:).';
    e1 = (p_tip - p_prev);
    e1 = e1 / norm(e1);             % local x1 (crack direction)
    e2 = [-e1(2); e1(1)];           % local x2 (normal to crack)

    % Global-from-local and local-from-global rotation matrices
    R_gl  = [e1, e2];               % x_global = R_gl * x_local
    R_loc = R_gl.';                 % x_local  = R_loc * x_global

    % ---------------------------------------------------------------------
    % 2.  BUILD HALF-CIRCLE IN LOCAL COORDINATES AND MAP TO GLOBAL
    % ---------------------------------------------------------------------
    nthet = 100;
    eps_th = 1e-3;                  % avoid crack line exactly

    theta = linspace(eps_th, pi-eps_th, nthet).';  % (0, π) upper half
    % Local coordinates of contour points (P: upper side)
    xP_loc = [ rI*cos(theta),  rI*sin(theta) ];    % [nthet x 2]
    xQ_loc = [ rI*cos(theta), -rI*sin(theta) ];    % mirror across x2=0

    % Map to global coordinates
    xP_gl = (R_gl * xP_loc.').';                  % [nthet x 2]
    xQ_gl = (R_gl * xQ_loc.').';                  % [nthet x 2]
    xP_gl = xP_gl + x_tip.';                      % shift to tip
    xQ_gl = xQ_gl + x_tip.';

    % (Optional) visualize contour
    % figure; plot(xP_gl(:,1),xP_gl(:,2),'b-',xQ_gl(:,1),xQ_gl(:,2),'r-'); axis equal;

    % ---------------------------------------------------------------------
    % 3.  FIND ELEMENTS FOR P, Q POINTS (T3 PARENT MESH)
    % ---------------------------------------------------------------------
    TR = triangulation(connect3, coord3);

    [elemP, baryP] = pointLocation(TR, xP_gl(:,1), xP_gl(:,2));
    [elemQ, baryQ] = pointLocation(TR, xQ_gl(:,1), xQ_gl(:,2));

    if any(isnan(elemP)) || any(isnan(elemQ))
        error('SIF_LEFM_CIRCLE: some contour points lie outside the mesh.');
    end

    % ---------------------------------------------------------------------
    % 4.  LOOP OVER HALF-CIRCLE, BUILD LOCAL FIELDS AT P and Q
    % ---------------------------------------------------------------------
    J1 = zeros(nthet,1);    % JI integrand at P(θ)
    J2 = zeros(nthet,1);    % JII integrand at P(θ)

    for k = 1:nthet

        % ---------- P: upper side ----------
        xiP  = baryP(k,1:2).';      % use first two barycentric coords
        eP   = elemP(k);            % element index
        nodesP = connect(eP,:);     % 6 T6 nodes
        XP   = coord(nodesP,:);     % [6 x 2] global coords of nodes

        % Get B, Det, and shape gradients dNdx at P
        [~, ~, dNdxP] = BN_local(xiP, XP);

        % Nodal displacements (6 x 2)
        uelP = [ U(2*nodesP-1), U(2*nodesP) ];

        % Global displacement gradient at P
        dux_dxP = dNdxP(1,:)*uelP(:,1);
        dux_dyP = dNdxP(2,:)*uelP(:,1);
        duy_dxP = dNdxP(1,:)*uelP(:,2);
        duy_dyP = dNdxP(2,:)*uelP(:,2);

        GradP_gl = [dux_dxP, dux_dyP;
                    duy_dxP, duy_dyP];      % [2 x 2]

        % Transform gradient to local crack coordinates
        GradP_loc = R_loc * GradP_gl * R_gl;     % [2 x 2]

        % Local strains at P
        eps11P = GradP_loc(1,1);
        eps22P = GradP_loc(2,2);
        gam12P = GradP_loc(1,2) + GradP_loc(2,1); % 2*eps12
        epsP   = [eps11P; eps22P; gam12P];

        % Local stresses at P (Dmat is isotropic plane strain)
        sigP = Dmat * epsP;    % [σ11; σ22; σ12]

        % Directional derivative ∂u^local/∂x1(local)
        du_dx1_P = GradP_loc(:,1);  % [du1/dx1; du2/dx1] at P

        % ---------- Q: lower side (mirror) ----------
        xiQ  = baryQ(k,1:2).';
        eQ   = elemQ(k);
        nodesQ = connect(eQ,:);
        XQ   = coord(nodesQ,:);

        [~, ~, dNdxQ] = BN_local(xiQ, XQ);

        uelQ = [ U(2*nodesQ-1), U(2*nodesQ) ];

        dux_dxQ = dNdxQ(1,:)*uelQ(:,1);
        dux_dyQ = dNdxQ(2,:)*uelQ(:,1);
        duy_dxQ = dNdxQ(1,:)*uelQ(:,2);
        duy_dyQ = dNdxQ(2,:)*uelQ(:,2);

        GradQ_gl = [dux_dxQ, dux_dyQ;
                    duy_dxQ, duy_dyQ];

        GradQ_loc = R_loc * GradQ_gl * R_gl;

        eps11Q = GradQ_loc(1,1);
        eps22Q = GradQ_loc(2,2);
        gam12Q = GradQ_loc(1,2) + GradQ_loc(2,1);
        epsQ   = [eps11Q; eps22Q; gam12Q];

        sigQ = Dmat * epsQ;

        du_dx1_Q = GradQ_loc(:,1);  % [du1/dx1; du2/dx1] at Q

        % -----------------------------------------------------------------
        % 5.  MODE SEPARATION IN LOCAL COORDINATES (Ishikawa–Kitagawa–Okamura)
        % -----------------------------------------------------------------
        % Strain vectors: [eps11; eps22; gamma12] where gamma12 = 2*eps12
        % Mode I:
        epsI = 0.5 * [ epsP(1) + epsQ(1);
                       epsP(2) + epsQ(2);
                       epsP(3) - epsQ(3) ];
        % Mode II:
        epsII = 0.5 * [ epsP(1) - epsQ(1);
                        epsP(2) - epsQ(2);
                        epsP(3) + epsQ(3) ];

        % Stress vectors:
        sigI = 0.5 * [ sigP(1) + sigQ(1);
                       sigP(2) + sigQ(2);
                       sigP(3) - sigQ(3) ];
        sigII = 0.5 * [ sigP(1) - sigQ(1);
                        sigP(2) - sigQ(2);
                        sigP(3) + sigQ(3) ];

        % Directional gradients ∂u/∂x1 (local) for modes I and II:
        % from u^I = (u + u')/2 (x1 component),
        %           (u - u')/2 (x2 component), etc.
        duI = 0.5 * [ du_dx1_P(1) + du_dx1_Q(1);   % du1^I/dx1
                      du_dx1_P(2) - du_dx1_Q(2) ]; % du2^I/dx1

        duII = 0.5 * [ du_dx1_P(1) - du_dx1_Q(1);  % du1^II/dx1
                       du_dx1_P(2) + du_dx1_Q(2) ];% du2^II/dx1

        % Energy densities
        UI  = 0.5 * (sigI.'  * epsI);
        UII = 0.5 * (sigII.' * epsII);

        % -----------------------------------------------------------------
        % 6.  J-INTEGRAND FOR JI AND JII AT POINT P (local coords)
        % -----------------------------------------------------------------
        % Local outward normal at P: n = [cosθ; sinθ]
        n1 = cos(theta(k));
        n2 = sin(theta(k));

        % Mode I traction at P: t = σ^I * n
        tI1 = sigI(1)*n1 + sigI(3)*n2;
        tI2 = sigI(3)*n1 + sigI(2)*n2;

        % JI integrand fI(θ) (in local coords, component along x1)
        J1(k) = UI*n1 - (tI1*duI(1) + tI2*duI(2));

        % Mode II traction at P
        tII1 = sigII(1)*n1 + sigII(3)*n2;
        tII2 = sigII(3)*n1 + sigII(2)*n2;

        % JII integrand fII(θ)
        J2(k) = UII*n1 - (tII1*duII(1) + tII2*duII(2));
    end

    % ---------------------------------------------------------------------
    % 7.  INTEGRATE ALONG HALF-CIRCLE AND MULTIPLY BY 2
    % ---------------------------------------------------------------------
    % ds = rI dθ, and we double for symmetric lower half.
    JI  = 2 * rI * trapz(theta, J1);
    JII = 2 * rI * trapz(theta, J2);

    % ---------------------------------------------------------------------
    % 8.  CONVERT TO KI, KII (plane strain)
    % ---------------------------------------------------------------------
    Eeff = E/(1-nu^2);   % plane strain modulus

    KI = sqrt(abs(JI) * Eeff);

    signK2 = 1;
    if JII < 0
        signK2 = -1;
    end
    KII = signK2 * sqrt(abs(JII) * Eeff);
end
