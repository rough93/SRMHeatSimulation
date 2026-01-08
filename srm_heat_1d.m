function srm_heat_1d()
clc;clear;

% user inputs
% materials
prop.name = 'APCP';     
prop.k = 0.35;  
prop.rho = 1700; 
prop.cp = 1200;

liner.name = 'EPDM';    
liner.k = 0.20; 
liner.rho = 1100; 
liner.cp = 1500;

caseM.name = 'Al6061';  
caseM.k = 167;  
caseM.rho = 2700; 
caseM.cp = 900;

caseCF.name = 'CF_composite';
caseCF.k = 1.0;
caseCF.rho = 1600;
caseCF.cp = 900;


% Layers & burn characteristics for sim
layers = struct([]);
layers(1).material = prop;   
layers(1).thickness = 0.0143;   % 14.3 mm propellant

layers(2).material = liner;  
layers(2).thickness = 0.002;   % 2 mm liner

layers(3).material = caseCF;  
layers(3).thickness = 0.00338;   % 3 mm case

geom.r_inner0 = 0.050;   % initial inner radius [m]
% burn.rburn = 1.0e-3;   % regression rate
burn.t_burn = 6; % layers(1).thickness / burn.rburn;      % burn duration, prescribed by test rather than calculated
burn.rburn = layers(1).thickness/burn.t_burn;

% Boundary conditions
BC.T_inner = 2300;       % combustion temperature kelvin
BC.T_amb = 300;        % ambient temp kelvin
BC.h_conv = 30;         % convection coefficient
BC.eps_rad = 0.8;        % emissivity
BC.sigma = 5.670374419e-8; % Stefan-Boltzmann
BC.T_lin = 800;        % linearization temperature for radiation [K]
BC.T_inner_off = BC.T_amb;  % drop inner surface temp after burnout to ambient
time.dt = 0.010;
time.t_final = 20.0;

% data output fidelity
cfg.store_every = 1;
cfg.max_nodes = 20000; % check number of nodes before run so my computer doesn't blow up.

% build 1-D grid
dx = burn.rburn * time.dt; % forcing dx to align with regression rate

[r_nodes, layer_idx_of_node] = build_grid_with_interfaces_on_nodes(geom, layers, dx);
n = numel(r_nodes);
r_nodes_initial = r_nodes;

if n > cfg.max_nodes
    error('Too many nodes (%d), increase dt.', n);
end

% give node calculation values
[k_node, rho_cp_node] = props_per_node(layers, layer_idx_of_node);

% assign h_effective (conv + linearized radiation)
Tlin = BC.T_lin;
h_eff = BC.h_conv + 4 * BC.eps_rad * BC.sigma * Tlin^3;

% initial condition
T = BC.T_amb * ones(n,1);

% Fo calculation to check later
alpha_prop = prop.k  / (prop.rho  * prop.cp);
Fo_prop = alpha_prop  * time.dt / dx^2;

alpha_liner = liner.k / (liner.rho * liner.cp);
Fo_liner = alpha_liner * time.dt / dx^2;

alpha_case = caseCF.k / (caseCF.rho * caseCF.cp);
Fo_case = alpha_case  * time.dt / dx^2;

fprintf('Fourier numbers (using dx = r_burn*dt): Fo_prop=%.3g, Fo_liner=%.3g, Fo_case=%.3g\n', Fo_prop, Fo_liner, Fo_case);

% time step loop
nsteps = round(time.t_final / time.dt);
fprintf('nodes=%d, dx=%.6g m, dt=%.6g s, steps=%d, h_eff=%.2f W/m^2-K\n', n, dx, time.dt, nsteps, h_eff);

times   = [];
Tsnap_r = {};
Tsnap_T = {};

for step = 0:nsteps
    t = step * time.dt;

    % store snapshot
    if mod(step, cfg.store_every) == 0
        times(end+1,1) = t;
        Tsnap_r{end+1} = r_nodes;
        Tsnap_T{end+1} = T;
    end

    if step == nsteps
        break;
    end

    % check if burning to set inner BC
    isBurning = (t <= burn.t_burn);
    if isBurning
        T_inner_curr = BC.T_inner;
    else
        T_inner_curr = BC.T_inner_off;
    end

    % assemble and solve
    [A, b] = assemble_BE_step(T, r_nodes, k_node, rho_cp_node, dx, time.dt, BC, h_eff, layer_idx_of_node, isBurning, T_inner_curr);
    T = A \ b;

    % advance burn front one node per step
    if (t < burn.t_burn) && (layer_idx_of_node(1) == 1)
        % remove inner node
        r_nodes(1)           = [];
        layer_idx_of_node(1) = [];
        T(1)                 = [];
        % recompute properties
        [k_node, rho_cp_node] = props_per_node(layers, layer_idx_of_node);
    end
end

% -------------------- plotting (2D sim + 1D-style at mid-z) --------------------
if ~exist('r','var') || ~exist('z','var')
    error('Plotting expects 2D variables r and z. Your run did not create them (still 1D/hybrid).');
end

Nt = numel(times);
Nr = numel(r);
Nz = numel(z);

iz_mid = round(Nz/2);
z_mid = z(iz_mid);

r_rel = r - geom.r_inner0;

% 1) Outer surface temperature vs time at z_mid
T_outer = zeros(Nt,1);
for n = 1:Nt
    Tn = Tsnap_T{n};            % Nr x Nz
    T_outer(n) = Tn(end, iz_mid);
end
figure;
plot(times, T_outer, 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('Outer Surface Temperature [K]');
title(sprintf('Case Outer Surface Temperature vs Time (z = %.3f m)', z_mid));

% 2) Radius vs time "3D" surface at z_mid with burned region filled for display
TplotMat = nan(Nt, Nr);
for n = 1:Nt
    t_n = times(n);
    r_front = geom.r_inner0 + burn.rburn * min(t_n, burn.t_burn);

    Tn = Tsnap_T{n};
    Tcol = Tn(:, iz_mid);

    burnedMask = (r <= r_front);
    unburnMask = ~burnedMask;

    if t_n <= burn.t_burn
        T_fill = BC.T_aw;     % or BC.T_inner if you still use it
    else
        T_fill = BC.T_amb;    % or BC.T_inner_off if desired
    end

    TplotMat(n, burnedMask) = T_fill;
    TplotMat(n, unburnMask) = Tcol(unburnMask);
end

[X, Y] = meshgrid(r_rel, times);

figure;
surf(X, Y, TplotMat, 'FaceColor','interp', 'EdgeColor','none');
colormap(turbo); colorbar;
grid on; box on;
xlabel('Radius from Original Inner Surface [m]');
ylabel('Time [s]');
zlabel('Temperature [K]');
title(sprintf('Temperature History at Mid-Axial Slice (z = %.3f m)', z_mid));

% red interface planes
hold on;
prop_thick  = layers(1).thickness;
liner_thick = layers(2).thickness;
interfaces = [prop_thick, prop_thick + liner_thick];

Tmin = min(TplotMat(:));
Tmax = max(TplotMat(:));
t_min = min(times);
t_max = max(times);

for kInt = 1:numel(interfaces)
    r_int = interfaces(kInt);
    Xp = r_int * ones(2,2);
    Yp = [t_min t_max; t_min t_max];
    Zp = [Tmin Tmin; Tmax Tmax];
    surf(Xp, Yp, Zp, 'FaceColor',[1 0 0], 'FaceAlpha',0.12, ...
         'EdgeColor',[1 0 0], 'EdgeAlpha',0.4, 'LineWidth',1.0);
end
hold off;

% 3) Contours (radius vs time) at z_mid
figure;
contourf(X, Y, TplotMat, 30, 'LineColor','none');
colormap(turbo); colorbar;
xlabel('Radius from Original Inner Surface [m]');
ylabel('Time [s]');
title(sprintf('Temperature Evolution (radius vs time) at z = %.3f m', z_mid));

% 4) Temperature vs time at selected radii at z_mid
r_sample_rel = [0.25, 0.50, 0.75] * r_rel(end);

figure; hold on;
plot(times, TplotMat(:, end), 'LineWidth', 2.0, 'DisplayName', 'Outer surface');
for kk = 1:numel(r_sample_rel)
    [~, idxr] = min(abs(r_rel - r_sample_rel(kk)));
    plot(times, TplotMat(:, idxr), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('r = %.3f m', r_rel(idxr)));
end
hold off; grid on;
xlabel('Time [s]');
ylabel('Temperature [K]');
title(sprintf('Temperature vs Time at Selected Radii (z = %.3f m)', z_mid));
legend('Location','best');

% 5) NEW: outer wall temperature over time along full axial length (z vs time)
Tout_z_t = nan(Nt, Nz);
for n = 1:Nt
    Tn = Tsnap_T{n};
    Tout_z_t(n,:) = Tn(end,:);
end

[Zgrid, Tgrid] = meshgrid(z, times);
figure;
surf(Zgrid, Tgrid, Tout_z_t, 'FaceColor','interp', 'EdgeColor','none');
colormap(turbo); colorbar;
grid on; box on;
xlabel('Axial Length z [m]');
ylabel('Time [s]');
zlabel('Outer Wall Temperature [K]');
title('Outer Wall Temperature vs Axial Length and Time');

end

%% -------------------- Helper functions --------------------

function [r_nodes, layer_idx_of_node] = build_grid_with_interfaces_on_nodes(geom, layers, dx)
% build radial node coordinates with layer interfaces laying on nodes.
r_inner = geom.r_inner0;
total_thickness = 0;
for L = layers
    total_thickness = total_thickness + L.thickness;
end

interfaces = zeros(numel(layers)+1,1);
interfaces(1) = r_inner;
acc = r_inner;
for i = 1:numel(layers)
    acc = acc + layers(i).thickness;
    interfaces(i+1) = acc;
end

r_nodes = r_inner;
layer_idx_of_node = 1;
for li = 1:numel(layers)
    rL = interfaces(li);
    rR = interfaces(li+1);
    r = rL;
    while true
        r_next = r + dx;
        if r_next >= rR - 1e-12
            if abs(rR - r_nodes(end)) > 1e-12
                r_nodes(end+1,1) = rR;
                layer_idx_of_node(end+1,1) = li;
            end
            break;
        else
            r_nodes(end+1,1) = r_next;
            layer_idx_of_node(end+1,1) = li;
            r = r_next;
        end
    end
end

end

function [k_node, rho_cp_node] = props_per_node(layers, layer_idx_of_node)
% assign node values based on location
n = numel(layer_idx_of_node);
k_node = zeros(n,1);
rho_cp_node = zeros(n,1);

for i = 1:n-1
    li = layer_idx_of_node(i);
    mat = layers(li).material;
    k_node(i) = mat.k;
    rho_cp_node(i) = mat.rho * mat.cp;
end

% set last node manually
mat = layers(end).material;
k_node(end) = mat.k;
rho_cp_node(end) = mat.rho * mat.cp;
end

function [A, b] = assemble_BE_step(Tprev, r_nodes, k_node, rho_cp_node, dx, dt, BC, h_eff, ~, isBurning, T_inner_curr)
% Assemble Backward Euler system A*T^{n+1} = b for one step
n = numel(r_nodes);
A = zeros(n,n);
b = zeros(n,1);

% Dirichlet inner boundary
if isBurning
    % prescribed temperature during burn
    A(1,1) = 1.0;
    b(1) = T_inner_curr;
else
    % q = 0 after burnout
    k_face = k_node(1);
    theta1 = rho_cp_node(1) * (dx/2) / dt;

    A(1,1) = theta1 + k_face/dx;
    A(1,2) = -k_face/dx;
    b(1) = theta1 * Tprev(1);
end

% interior nodes, backward euler equations
for i = 2:n-1
    % effective k to faces
    kL = 0.5*(k_node(i-1) + k_node(i)); % conductivity left
    kR = 0.5*(k_node(i) + k_node(i+1)); % conductivity right
    C_L = kL / dx^2;    % left conduction
    C_R = kR / dx^2;    % right conduction

    theta = rho_cp_node(i) * dx/dt;  % transient coefficient

    A(i,i-1) = -C_L;
    A(i,i) = theta + C_L + C_R;
    A(i,i+1) = -C_R;
    b(i) = theta * Tprev(i);
end

% Robin outside boundary
i = n;
k_face = k_node(i);
thetaN = rho_cp_node(i) * (dx/2) / dt;
A(i,i-1) = -k_face / dx;
A(i,i)   = thetaN + k_face/dx + h_eff;
b(i)     = thetaN * Tprev(i) + h_eff * BC.T_amb;

end
