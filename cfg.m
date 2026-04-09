function C = cfg()
% All parameters; nelvert is locked to 6

cm = .01;
C.A = 10*cm;             % width
C.B = 10*cm;             % height
C.a =  2*cm;             % initial crack length

C.nelvert = 6;           % enforced
C.fact = 20;             % visualization scale
C.eps1 = 1e-8;
C.hmax = C.B/40;
C.ncoh = 40;
C.hgrad = 1.1;
C.da    = 0.06*C.a;
C.chw = 1/8 * C.da / C.ncoh;      % channel width

% Single theta (deg), as in your snippet:
C.theta_deg = 20; %1.655; 3.390
C.alf1 = -20.0*pi/180;     % first leg
C.alf2 = C.theta_deg*pi/180;

% Local frame (rotation) at the last leg
calf = cos(C.alf2); salf = sin(C.alf2);
C.malf = [calf, salf; -salf, calf];

C.V0 = [0, 0];      % mouth location (global)
C.V1 = C.a*[cos(C.alf1), sin(C.alf1)];
C.V2 = C.V1 + C.da*[cos(C.alf2), sin(C.alf2)];

% Material
C.nu = 0.3; C.E2 = 4e3;
C.Dmat = C.E2*[[1-C.nu, C.nu, 0];
               [C.nu, 1-C.nu, 0];
               [0, 0, (1-2*C.nu)/2]]/(1+C.nu)/(1-2*C.nu);
C.G12 = C.E2/2/(1+C.nu);

% TSL params
C.phi  = [800, 1200];
%C.rphi = C.phi(1)/C.phi(2);
C.rphi = 0.4;
C.eta_BK = [1, 1.5, 3];
C.a1   = 0.002;
C.a2   = [.9, .5];
C.rEs  = 100;
C.sigmax = [C.E2, C.G12]/C.rEs;
c_area = (3 - 2*C.a1 + 3*C.a2)/6;             % vectorized
C.Dmax  = C.phi ./ c_area ./ (C.sigmax*1e6);


% load guess
C.sig_guess = C.sigmax(1)/8;

% Plot
C.plotGeom = false;
C.plotMesh = false;

end
