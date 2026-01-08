clc; clear;

% -------------------- Inputs --------------------
% Materials
mat.prop  = struct('name','APCP',        'k',0.35, 'rho',1700, 'cp',1200);
mat.liner = struct('name','EPDM',        'k',0.20, 'rho',1100, 'cp',1500);
mat.case  = struct('name','CF_composite','k',1.0,  'rho',1600, 'cp',900);

% Layers (radial only)
layers = struct([]);

layers(1).material  = mat.prop;
layers(1).thickness = 0.0143;

layers(2).material  = mat.liner;
layers(2).thickness = 0.0020;

layers(3).material  = mat.case;
layers(3).thickness = 0.00338;

% Geometry
geom = struct('r_inner0',0.050, 'Lz',0.10);

% Burn
burnCfg = struct();
burnCfg.t_burn = 6.0;
burnCfg.rburn  = layers(1).thickness / burnCfg.t_burn;

% Boundary conditions
BC = struct();
BC.T_amb   = 300;
BC.h_conv  = 30;
BC.eps_rad = 0.8;
BC.sigma   = 5.670374419e-8;
BC.T_lin   = 800;

% Inner gas-side Robin (during burn)
BC.h_gas = 1500;     % placeholder
BC.T_aw  = 2800;     % placeholder

% Optional cavity convection after burn
BC.h_cav = 20;       % example
BC.T_cav = BC.T_amb;

% Derived outer h_eff
BC.h_eff = BC.h_conv + 4 * BC.eps_rad * BC.sigma * BC.T_lin^3;

% Discretization
num = struct('Nr',250, 'Nz',100);

% Time
time = struct('dt',0.05, 't_final',40.0);

% Output / debug
cfg = struct();
cfg.store_every = 1;
cfg.max_nodes = 250000;
cfg.axial_adiabatic = true;

% -------------------- Build grid and material fields --------------------
% build_grid_rz should return grid struct with: r,z,R,Z,dr,dz,t_wall,r_outer
grid = build_grid_rz(geom, layers, num);

% material_map_r should return fields: k(Nr,Nz), rhocp(Nr,Nz), layerId(Nr,Nz),
% and optionally interfaces for plotting
matField = material_map_r(grid, mat, layers, geom);

% Diagnostics
if num.Nr * num.Nz > cfg.max_nodes
    warning('Grid has %d nodes; consider raising cfg.max_nodes or reducing Nr/Nz.', num.Nr*num.Nz);
end

nsteps = round(time.t_final / time.dt);
fprintf('Nr=%d Nz=%d nodes=%d | dr=%.3g dz=%.3g dt=%.3g steps=%d | h_eff=%.2f\n', ...
    num.Nr, num.Nz, num.Nr*num.Nz, grid.dr, grid.dz, time.dt, nsteps, BC.h_eff);

% Optional Fo diagnostics (use mid material properties or per-layer like before)
alpha_prop  = mat.prop.k   / (mat.prop.rho   * mat.prop.cp);
alpha_liner = mat.liner.k  / (mat.liner.rho  * mat.liner.cp);
alpha_case  = mat.case.k   / (mat.case.rho   * mat.case.cp);

Fo_r_prop  = alpha_prop  * time.dt / grid.dr^2;
Fo_r_liner = alpha_liner * time.dt / grid.dr^2;
Fo_r_case  = alpha_case  * time.dt / grid.dr^2;

Fo_z_prop  = alpha_prop  * time.dt / grid.dz^2;
Fo_z_liner = alpha_liner * time.dt / grid.dz^2;
Fo_z_case  = alpha_case  * time.dt / grid.dz^2;

fprintf('Fo_r: prop=%.3g liner=%.3g case=%.3g | Fo_z: prop=%.3g liner=%.3g case=%.3g\n', ...
    Fo_r_prop, Fo_r_liner, Fo_r_case, Fo_z_prop, Fo_z_liner, Fo_z_case);

% -------------------- Initialize state --------------------
T = BC.T_amb * ones(num.Nr, num.Nz);

% Snapshot storage
Tsnap_T = {};
times   = [];

% Mid-slice index used in plots
iz_mid = round(num.Nz/2);

% -------------------- Cached system --------------------
iFront_prev = -1;
A_prev  = [];
dA_prev = [];           % LU cache
BCinner_prev  = struct('h', nan, 'Tinf', nan);

% -------------------- Time loop --------------------
% -------------------- Time loop --------------------
iFront_prev    = -1;
BCinner_prev   = struct('h', nan, 'Tinf', nan);
A_prev         = [];
dA_prev        = [];

for step = 0:nsteps
    t = step * time.dt;

    % Burn front info (must provide burn.iFront, burn.isBurning, burn.isInnerFace)
    burn = burn_front_rzt(grid, burnCfg, t);

    % Inner BC switching (based on burn state)
    if burn.isBurning
        BCinner = struct('h', BC.h_gas, 'Tinf', BC.T_aw);
    else
        % post-burn: set BC.h_cav = 0 for adiabatic if desired
        BCinner = struct('h', BC.h_cav, 'Tinf', BC.T_cav);
    end

    % Rebuild when interface index changes or inner BC changes
    bcChanged   = (BCinner.h ~= BCinner_prev.h) || (BCinner.Tinf ~= BCinner_prev.Tinf);
    needRebuild = (burn.iFront ~= iFront_prev) || bcChanged || isempty(dA_prev);

    if needRebuild
        [A_prev, b] = assemble_BE_2d_axisym(T, grid, matField, BC, burn, BCinner, time.dt, cfg);
        dA_prev = decomposition(A_prev, 'lu');
        iFront_prev  = burn.iFront;
        BCinner_prev = BCinner;
    else
        b = assemble_b_only_axi2d(T, grid, matField, BC, burn, BCinner, time.dt, cfg);
    end

    % Solve
    Tvec = dA_prev \ b;

    % Scatter back into full grid (cavity filled for convenience)
    Tnew = BC.T_amb * ones(num.Nr, num.Nz);
    Tnew(burn.isSolid) = Tvec;
    T = Tnew;

    % ---- Snapshot storage (THIS WAS MISSING) ----
    if step == 0 || mod(step, cfg.store_every) == 0 || step == nsteps
        times(end+1,1)     = t;
        Tsnap_T{end+1,1}   = T;
    end

    if mod(step,50)==0
        fprintf('step=%d/%d t=%.2f r_front=%.5f iFront=%d | Tmin=%.2f Tmax=%.2f TouterMid=%.2f\n', ...
            step, nsteps, t, burn.r_front, burn.iFront, min(T(:)), max(T(:)), T(end, round(num.Nz/2)));
    end
end


fprintf('DEBUG: nsteps=%d, stored Nt=%d, Tsnap_T cells=%d\n', nsteps, numel(times), numel(Tsnap_T));

% -------------------- Pack sim struct for debugging/replotting --------------------
sim = struct();
sim.geom = geom;
sim.layers = layers;
sim.mat = mat;
sim.matField = matField;
sim.burnCfg = burnCfg;
sim.BC = BC;
sim.num = num;
sim.time = time;
sim.cfg = cfg;
sim.grid = grid;

sim.out = struct();
sim.out.times = times;
sim.out.Tsnap_T = Tsnap_T;
sim.out.iz_mid = iz_mid;

% -------------------- Plot --------------------
postprocess_plots_2d(sim);

