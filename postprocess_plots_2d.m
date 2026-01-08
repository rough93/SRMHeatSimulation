function postprocess_plots_2d(sim)
%POSTPROCESS_PLOTS_2D 1D-style plots at mid-z + outer wall z–t plot
% Assumes snapshots sim.out.Tsnap_T are Nr x Nz each.

g      = sim.grid;
geom   = sim.geom;
layers = sim.layers;
BC     = sim.BC;
burnCfg = sim.burnCfg;

times   = sim.out.times(:);      % FORCE Nt x 1
Tsnap_T = sim.out.Tsnap_T;

Nt = numel(times);
Nr = numel(g.r);
Nz = numel(g.z);

% Mid-z index
if isfield(sim.out,'iz_mid') && ~isempty(sim.out.iz_mid)
    iz_mid = sim.out.iz_mid;
else
    iz_mid = round(Nz/2);
end
iz_mid = max(1, min(Nz, iz_mid));
z_mid  = g.z(iz_mid);

% Force r_rel to be 1 x Nr for meshgrid consistency
r_rel = (g.r - geom.r_inner0).';  % 1 x Nr

% -------- preallocate derived arrays --------
T_outer_mid = zeros(Nt,1);        % Nt x 1
T_outer_z_t = nan(Nt, Nz);        % Nt x Nz
T_mid_r_t   = nan(Nt, Nr);        % Nt x Nr

% -------- build derived arrays from snapshots --------
for n = 1:Nt
    Tn = Tsnap_T{n};
    if ~isequal(size(Tn), [Nr, Nz])
        error('Snapshot %d size is %dx%d but expected %dx%d.', ...
            n, size(Tn,1), size(Tn,2), Nr, Nz);
    end

    T_outer_mid(n)   = Tn(end, iz_mid);
    T_outer_z_t(n,:) = Tn(end, :);
    T_mid_r_t(n,:)   = Tn(:, iz_mid).';   % row 1xNr
end

% ---------------- 1) Outer surface temperature vs time at mid-z ----------------
figure;
plot(times, T_outer_mid, 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('Outer Surface Temperature [K]');
title(sprintf('Case Outer Surface Temperature vs Time (z = %.3f m)', z_mid));
set(gca,'YDir','normal');

% ---------------- 2) Radius vs time surface at mid-z (with burned region filled) ----------------
TplotMat = nan(Nt, Nr);  % Nt x Nr

for n = 1:Nt
    t_n = times(n);
    r_front = geom.r_inner0 + burnCfg.rburn * min(t_n, burnCfg.t_burn);

    Trow = T_mid_r_t(n,:);  % 1 x Nr

    burnedMask = (g.r <= r_front).';   % 1 x Nr
    unburnMask = ~burnedMask;

    if t_n <= burnCfg.t_burn
        T_fill = BC.T_aw;
    else
        T_fill = BC.T_amb;
    end

    TplotMat(n, burnedMask) = T_fill;
    TplotMat(n, unburnMask) = Trow(unburnMask);
end

% Ensure surf inputs are exactly Nt x Nr
TplotMat = reshape(TplotMat, [Nt, Nr]);

% Build matching grids (Nt x Nr)
[X, Y] = meshgrid(r_rel, times);

% Final assert (optional but helpful)
if ~isequal(size(X), size(TplotMat)) || ~isequal(size(Y), size(TplotMat))
    error('surf mismatch: size(X)=%dx%d size(Y)=%dx%d size(TplotMat)=%dx%d', ...
        size(X,1), size(X,2), size(Y,1), size(Y,2), size(TplotMat,1), size(TplotMat,2));
end

[X, Y] = meshgrid(r_rel, times);   % X,Y are Nt x Nr

figure;
surf(X, Y, TplotMat, 'FaceColor','interp', 'EdgeColor','none');
colormap(turbo); colorbar;
grid on; box on;
xlabel('Radius from Original Inner Surface [m]');
ylabel('Time [s]');
zlabel('Temperature [K]');
title(sprintf('Temperature History at Mid-Axial Slice (z = %.3f m)', z_mid));
set(gca,'YDir','normal');
set(gca,'XDir','reverse');
view(45,30);

% Interface planes in r_rel coordinates
hold on;
prop_thick  = layers(1).thickness;
liner_thick = layers(2).thickness;
interfaces_rel = [prop_thick, prop_thick + liner_thick];

Tmin = min(TplotMat(:));
Tmax = max(TplotMat(:));
t_min = min(times);
t_max = max(times);

for kInt = 1:numel(interfaces_rel)
    x_int = interfaces_rel(kInt);
    surf(x_int*ones(2,2), [t_min t_max; t_min t_max], [Tmin Tmin; Tmax Tmax], ...
        'FaceColor',[1 0 0], 'FaceAlpha',0.12, ...
        'EdgeColor',[1 0 0], 'EdgeAlpha',0.4, 'LineWidth',1.0);
end
hold off;

% ---------------- 3) Contour plot (radius vs time) ----------------
figure;
contourf(X, Y, TplotMat, 30, 'LineColor','none');
colormap(turbo); colorbar;
xlabel('Radius from Original Inner Surface [m]');
ylabel('Time [s]');
title(sprintf('Temperature Evolution (radius vs time) at z = %.3f m', z_mid));

% ---------------- 4) Temperature vs time at selected radii ----------------
r_queries_rel = [0.25, 0.50, 0.75] * (g.r(end) - geom.r_inner0);
r_queries_abs = geom.r_inner0 + r_queries_rel;

figure; hold on;
plot(times, T_outer_mid, 'LineWidth', 2.0, 'DisplayName', 'Outer surface');

for kk = 1:numel(r_queries_abs)
    T_here = sample_radius_timeseries(Tsnap_T, times, g, geom, burnCfg, iz_mid, r_queries_abs(kk));
    plot(times, T_here, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('r_{rel} = %.4f m', r_queries_rel(kk)));
end

hold off; grid on;
xlabel('Time [s]'); ylabel('Temperature [K]');
title(sprintf('Temperature vs Time at Selected Radii (z = %.3f m)', z_mid));
legend('Location','best');

% ---------------- 5) Outer wall temperature vs axial length and time ----------------
[Zgrid, Tgrid] = meshgrid(g.z, times);  % Nt x Nz
figure;
surf(Zgrid, Tgrid, T_outer_z_t, 'FaceColor','interp', 'EdgeColor','none');
colormap(turbo); colorbar;
grid on; box on;
xlabel('Axial Length z [m]');
ylabel('Time [s]');
zlabel('Outer Wall Temperature [K]');
title('Outer Wall Temperature vs Axial Length and Time');
set(gca,'YDir','normal');
view(45,30);

end

function Tseries = sample_radius_timeseries(Tsnap_T, times, g, geom, burnCfg, iz, r_query)
Nt = numel(times);
Tseries = nan(Nt,1);

for n = 1:Nt
    t = times(n);
    r_front = geom.r_inner0 + burnCfg.rburn * min(t, burnCfg.t_burn);

    if r_query <= r_front
        Tseries(n) = NaN;  % burned-out region
    else
        [~, ir] = min(abs(g.r - r_query));
        Tseries(n) = Tsnap_T{n}(ir, iz);
    end
end
end
