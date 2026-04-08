function plot_stress_subfigure(ax, mesh, U, cz, mat, quad, nstress, fact)
%PLOT_STRESS_SUBFIGURE Plot stress field into axes ax.
%
% Based on the first part of DispStressExt(...), but without creating a new figure.
%
% Inputs:
%   ax      - target axes
%   mesh,U,cz,mat,quad,nstress,fact
%             same meaning as in DispStressExt

% --- required inputs ---
Dmat = mat.D;
E    = mat.E;
nu   = mat.nu;

% --- Stress field ---
[coord1, sig_node] = StressExt(mesh, U, mat, quad, fact);
coord   = mesh.coord;
connect = mesh.connect;

xlim([-.02 0.32]);
ylim([-.24 0.24]);

nelem = size(connect,1);
vert  = [1,4,2,5,3,6];   % node order for plotting T6

% normalization for color field
if nstress == 3
    smax = max(1e-12, cz.sigmax(2));   % shear component normalization
elseif nstress == 4
    smax = 1;                          % von Mises absolute
else
    smax = max(1e-12, cz.sigmax(1));   % normal component normalization
end

X = zeros(nelem, numel(vert));
Y = X;
Z = X;

for e = 1:nelem
    sig0 = sig_node(connect(e,vert), :);   % [sxx, syy, sxy]

    if nstress < 4
        Z(e,:) = sig0(:,nstress) / smax;
    else
        % 3D-equivalent von Mises (plane strain sigma_z)
        eps0 = (Dmat \ sig0')';
        sigx = sig0(:,1);
        sigy = sig0(:,2);
        sigxy = sig0(:,3);
        sigz = E*nu/(1+nu)/(1-2*nu) * (eps0(:,1) + eps0(:,2));
        vm   = sqrt((sigx - sigy).^2 + (sigy - sigz).^2 + (sigz - sigx).^2 + 6*sigxy.^2)/sqrt(2);
        Z(e,:) = vm;
    end

    X(e,:) = coord1(connect(e,vert),1);
    Y(e,:) = coord1(connect(e,vert),2);
end

% --- plotting ---
cla(ax);
hold(ax,'on');
view(ax,2);
set(ax,'DataAspectRatio',[1 1 1]);

col = linspace(0,1,21)';
cmap = jet(numel(col)-1);

if nstress < 4
    clim(ax, [col(1) col(end)]);
end
colormap(ax, cmap);

patch(ax, X', Y', Z', ...
    'EdgeColor', 'none', ... %0.6*[1,1,1]
    'LineWidth', 0.05, ...
    'FaceColor', 'interp');

cb = colorbar(ax);
if nstress < 4
    cb.Ticks = col;
end

if nstress == 1
    title(ax, '$\sigma_{xx} / \sigma_{(n)}^{\max}$', 'Interpreter','latex');
elseif nstress == 2
    title(ax, '$\sigma_{yy} / \sigma_{(n)}^{\max}$', 'Interpreter','latex');
elseif nstress == 3
    title(ax, '$\sigma_{xy} / \sigma_{(t)}^{\max}$', 'Interpreter','latex');
else
    title(ax, 'von Mises (absolute)');
end

box(ax,'on');
end