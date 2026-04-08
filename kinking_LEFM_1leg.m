function [theta_k, theta_deg, KI, KII] = kinking_LEFM_1leg(C,sigma0)
%KINKING_LEFM_1LEG  Compute LEFM kinking angle for a 1-leg crack.
%
%   [theta_k, theta_deg, KI, KII] = kinking_LEFM_1leg(C)
%
%   Required fields of C (adapt to your notation):
%       A, B      - strip half-widths: x in [0,A], y in [-B,B]
%       a         - crack length (mouth -> tip)
%       alf1      - crack angle (radians) from +x axis
%       chw       - half-width of refined "pencil" channel around crack
%       hmax      - target max element size away from crack
%       hgrad     - mesh gradation
%       eps1      - smoothing for edge loads (as in your code)
%       E2, nu    - elastic constants (plane strain)
%       G12       - shear modulus (if needed by your Dmat builder)
%       Dmat      - 3x3 elasticity matrix for plane strain
%
%   This solves the purely elastic problem with unit nominal stress
%   on the top/bottom edges (sigma = 1 MPa) and extracts KI, KII at
%   the tip of the straight crack, then computes the MTS kinking angle.

    % ---------- geometry + mesh + stiffness ----------
    [ mesh, mat, quad, elod, Stif, fixvar, id_tip] = geom_pencil_1leg(C,0.005);

    coord   = mesh.coord;
    connect = mesh.connect;

    nnod = size(coord,1);

    F = zeros(2*nnod,1);           % global load vector

    if ~isempty(elod)
        ids = elod(:,1);           % node indices
        w   = elod(:,2);           % 1D shape-integrated weights along edges

        % vertical nodal forces: Fy = sigma0 * w
        Fy = sigma0 * w;

        % x-component is zero; add Fy to DOF 2*id
        F(2*ids) = F(2*ids) + Fy;
    end


    % Apply homogeneous Dirichlet BCs on fixed DOFs
    F(fixvar) = 0;

    % ---------- solve linear system ----------
    U = Stif \ F;   % displacement vector (2*nnod x 1)

    if C.plotMesh
        cz = struct('sigmax',sigma0,'a',C.a,'tip_id',id_tip);
        DispStressLEFM(mesh, U, cz, mat, quad, 2, 20);
    end

    % ---------- SIFs at the tip ----------
    rI     = 0.5*C.a;                 % J-contour radius
    mat = struct('E',C.E2,'nu',C.nu,'Dmat',C.Dmat);
    V = [C.V0; C.V1];

    [KI, KII] = SIF_LEFM_circle2(mesh, U, V, mat, rI);

    % ---------- MTS kinking angle ----------
    [theta_k, theta_deg] = kink_angle_LEFM_MTS(KI, KII);

    fprintf('LEFM SIFs: KI = %.4g, KII = %.4g (units MPa*sqrt(m))\n', KI, KII);
    fprintf('MTS kink angle: %.3f rad = %.2f deg\n', theta_k, theta_deg);
end
