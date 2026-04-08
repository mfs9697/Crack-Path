%% main_crit3.m  -- driver for critical state solve using f1_crit3 (no tip cohesive elem)

clear; clc;

%% ---------- Stage 1: config + geometry ----------
C = cfg_min_geom();   % geometry + CZ params + solver params

% Geometry & mesh (+materials/quad, used for post)
[mesh, mat, quad, elod, Stif, cohes] = geom_pencil(C);

% Cohesive-zone bundle (keep same fields as before)
cz = struct('cohes',cohes, 'Q',C.Q,...
    'Dmax',C.Dmax, 'sigmax',C.sigmax, ...
    'a1',C.a1, 'a2',C.a2, 'rphi',C.rphi);

% Optional handle (only for diagnostics/plots; f1_crit3 calls tsl2 directly)
cz.tsl = @(DnDt) tsl2(DnDt, cz);

% Plotting / post flags
opts = struct('fact', C.fact, 'nstress', 2);

%% ---------- Stage 2: solve augmented system [u; sig] ----------
ndof = 2*size(mesh.coord,1);

% Jacobian sparsity (bulk + CZ band + 1 col/row for sigma)
Jpat = spones(Stif);

% Cohesive station list to use (prefer 'cohes' argument passed to f1_crit3)
% cohes is (2*ncoh+1) x 2  [upper lower], interleaved v,m,v..., mouth->tip
mcoh      = size(cohes,1);
ncoh_elem = (mcoh-1)/2;    % total cohesive quadratic elements

% Add CZ coupling pattern for ALL cohesive elements (including the tip element)
for q = 1:ncoh_elem
    ind1 = 2*(q-1)+(1:3);                 % [vertex, mid, vertex] station rows
    up   = cohes(ind1,1);
    lo   = cohes(ind1,2);

    dof = reshape([2*up-1; 2*up; 2*lo-1; 2*lo], [], 1);
    Jpat(dof,dof) = 1;
end

% Augmented pattern (append sigma)
Jpat_aug = blkdiag(Jpat, speye(1));

% TypicalX scaling (keep your pragmatic scaling)
x_scale = max(1e-9, (C.sig_guess/max(cz.sigmax(1),1))*C.B) * ones(ndof+1,1);

opts_solve = optimoptions('fsolve', ...
    'Algorithm','trust-region-dogleg', ...
    'SpecifyObjectiveGradient', true, ...
    'JacobPattern', Jpat_aug, ...
    'ScaleProblem','jacobian', ...
    'TypicalX', x_scale, ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-12, ...
    'MaxIterations',100, ...
    'Display','iter-detailed');

% Initial guess
x = [zeros(ndof,1); C.sig_guess];

% Residual/Jacobian handle (any h)
fun = @(xx) f1_crit3(xx, mesh, cz, Stif, elod, cohes);

% The full nonlinearity 
[x, fval, flag, out] = fsolve(@(xx) fun(xx), x, opts_solve);
fprintf('||F||=%.3e, iters=%d, flag=%d\n', norm(fval), out.iterations, flag);

% Unpack
u_sol  = x(1:ndof);
sig_cr = x(end);
fprintf('sig_crit = %.10g\n', sig_cr);

%% ---------- Stage 3: post ----------
DispStressExt(mesh, u_sol, cz, mat, quad, opts.nstress, opts.fact);
