function SRM_Heat_3D_Polar_endcap()
clc;

% 2D transient heat conduction in polar cross-section (r,theta) for SRM wall
% Finite-volume, Backward Euler, moving inner boundary (burnback).
%
% Energy-consistent geometry:
%   CV volume per unit length L=1: V = 0.5*(r_{i+1/2}^2 - r_{i-1/2}^2)*dtheta
%   Radial face area per unit length: A_r = r_face * dtheta
%   Theta face area per unit length: A_theta = dr_i
%   Theta face distance: ds = r_i * dtheta

%% -------------------- User inputs (same intent as your 1D code) --------------------
% materials
prop.name = 'APCP';
prop.k = 0.35;
prop.rho = 1700;
prop.cp = 1200;
prop.burnrate = 0.002383; % [m/s]

liner.name = 'EPDM';
liner.k = 0.20;
liner.rho = 1100;
liner.cp = 1500;

caseCF.name = 'CF_composite';
caseCF.k = 1.0;
caseCF.rho = 1600;
caseCF.cp = 900;

caseAl.name = 'Aluminum';
caseAl.k   = 205;     % W/m-K
caseAl.rho = 2700;    % kg/m^3
caseAl.cp  = 900;     % J/kg-K

% layers
layers = struct([]);
layers(1).material = prop;
layers(1).thickness = 0.0824;
layers(2).material = liner;
layers(2).thickness = 0.004;
layers(3).material = caseAl;
layers(3).thickness = 0.0127;

% geometry / burn
geom.r_inner0 = 0.015;
burn.rburn = prop.burnrate;
burn.t_burn = layers(1).thickness / burn.rburn;   % [s]


cfg.store_every = 1;     % store every N steps
cfg.max_cells = 20e6;      % safety
cfg.Ntheta = 16;          % circumferential resolution (>= 8)
cfg.Nz = 600;              % axial resolution (>=1)
cfg.endBC = 'adiabatic';    % 'adiabatic' or 'robin' at z=0 and z=L

time.dt = 0.1;
time.t_final = 200.0;

% boundary conditions
BC.T_gas = 2300;          % "inner" gas temperature
BC.T_amb = 300;           % ambient
BC.h_out = 30;            % outer convection coefficient
BC.h_in  = 2000;          % inner convection coefficient (set realistically later)
BC.eps_in  = 0.8;
BC.eps_out = 0.8;
BC.sigma = 5.670374419e-8;

% ---- Step 3: aft cavity inner-wall gas-side BC (reduced vs port) ----
BC.h_in_cav   = 200;     % [W/m^2-K] cavity-side convection (lower than BC.h_in)
BC.eps_in_cav = 0.8;     % [-] emissivity for cavity region (can keep same as eps_in)
BC.T_gas_cav_mode = 'same';  % 'same' uses T_gas, or 'ambient' uses BC.T_amb

BC.h_z0 = 30;               % end convection at z=0 if endBC='robin'
BC.h_zL = 30;               % end convection at z=L if endBC='robin'
cfg.L = 1.66;              % unit length (m)



% ---- Step 4: nozzle/throat proxy for aft-cavity convection ----
cfg.noz.enableProxy = true;
cfg.noz.r_t = 0.010;                 % [m] throat radius (user input)
cfg.noz.CdA = 1.0;                   % dimensionless scale factor (leave 1 for now)
cfg.noz.h_cav_base = 200;            % [W/m^2-K] baseline cavity h (what you used in Step 3)
cfg.noz.alpha = 0.5;                 % exponent for area-ratio scaling (0.5 is a safe start)
cfg.noz.h_cav_min = 10;              % [W/m^2-K] floor to avoid "zero convection"
cfg.noz.h_cav_max = 800;             % [W/m^2-K] ceiling to avoid runaway h

% ---- Optional endcap / bulkhead thermal mass (lumped nodes coupled to z-ends) ----
cfg.endcap.enable = true;        % set false to disable endcap nodes
cfg.endcap.ends   = 'both';      % 'z0', 'zL', or 'both'
cfg.endcap.thickness = 0.012;    % [m] effective bulkhead thickness participating thermally
cfg.endcap.k      = 205;         % [W/m-K] bulkhead conductivity (e.g., aluminum ~205)
cfg.endcap.rho    = 2700;        % [kg/m^3]
cfg.endcap.cp     = 900;         % [J/kg-K]
cfg.endcap.h_in   = 0;           % [W/m^2-K]
cfg.endcap.eps_in = 0;           % [-]
cfg.endcap.h_out  = 0;           % [W/m^2-K]
cfg.endcap.eps_out = 0;          % [-]

% ---- Endcap gas exposure flags (Step 1: inhibited grain end face) ----
cfg.endcap.gas_exposed_z0 = false;   % inhibited/insulated end at z=0
cfg.endcap.gas_exposed_zL = true;    % nozzle/free-volume end at z=L

% ---- Step 2: axial insulation on inner wall near inhibited end ----
cfg.wall_insul.enable = true;
cfg.wall_insul.end    = 'z0';      % 'z0' or 'zL'
cfg.wall_insul.length = 0.03;      % [m] axial insulation length
cfg.wall_insul.factor = 0.0;       % 0 = fully insulated, 0.1 = partial, etc.

% ---- Axial regression / free-volume model (minimal) ----
cfg.grain.enable_axial_regression = true;
cfg.grain.Lg0 = cfg.L;         % initial grain length, can be < cfg.L later
cfg.grain.rb_end = burn.rburn; % rb_end = rb_radial

% ---- runtime status monitoring ----
cfg.status_every = 20;     % print every N time steps
cfg.time_warn    = 2.0;    % warn if a step takes > this many seconds
cfg.doEnergyDiagnostics = false;     % master switch for all plots/prints below
cfg.doDeletionDiagnostics = false;   % requires Qcond12A_hist + willDel_hist

% ---- NEW: diag flag (used by the drop-in block) ----
cfg.diag.consistency = true;

%% -------------------- Build polar FV grid --------------------
dr = burn.rburn * time.dt;         % radial step aligned to burn per step
[dtheta, theta_centers] = uniform_theta_grid(cfg.Ntheta);

[r_faces, r_centers, layer_idx_of_ring] = build_radial_faces_with_interfaces(geom, layers, dr);
Nr = numel(r_centers);
Ntheta = cfg.Ntheta;

Nz = cfg.Nz;
dz = cfg.L / Nz;

nCells = Nr * Ntheta * Nz;  % 3D cells
if nCells > cfg.max_cells
    error('Too many cells (%d). Reduce Ntheta or increase dt.', nCells);
end

% material props per radial ring (piecewise constant in theta by default)
[k_ring, rhoCp_ring] = ring_props(layers, layer_idx_of_ring);

% initial condition
T = BC.T_amb * ones(Nr, Ntheta, Nz);

% Endcap lumped temperatures (if enabled)
if isfield(cfg,'endcap') && isfield(cfg.endcap,'enable') && cfg.endcap.enable
    Tcap0 = BC.T_amb;
    TcapL = BC.T_amb;
else
    Tcap0 = [];
    TcapL = [];
end

% diagnostics: rough Fourier number using smallest ring dr
alpha_prop = prop.k / (prop.rho*prop.cp);
Fo = alpha_prop * time.dt / dr^2;
fprintf('Grid: Nr=%d, Ntheta=%d, dr=%.6g m, dt=%.6g s, Fo_prop~%.3g\n', Nr, Ntheta, dr, time.dt, Fo);

%% -------------------- Time integration --------------------
nsteps = round(time.t_final / time.dt);
timing = struct( ...
    'cache', 0.0, ...
    'assembly', 0.0, ...
    'prec', 0.0, ...
    'solve', 0.0, ...
    'diagAb', 0.0, ...
    'bcHeat', 0.0, ...
    'deletion', 0.0, ...
    'store', 0.0, ...
    'other', 0.0 );

timing_step = timing;   % optional: for per-step debugging

% ---- Energy balance tracking ----
N = nsteps;  % number of steps with solves
E_solid_hist   = nan(N+2,1);   % carry state, needs step+2 indexing
dE_hist        = nan(N,1);
dQnet_hist     = nan(N,1);
closure_hist   = nan(N,1);
relerr_hist    = nan(N,1);
Qcond12A_hist  = nan(N,1);
Erem_hist      = nan(N,1);
willDel_hist   = false(N,1);

closure_pre_hist  = nan(N,1);
closure_post_hist = nan(N,1);
Qdel_hist         = nan(N,1);

% ---- Ab-implied + residual storage ----
dQimplied_hist      = nan(N,1);
dQimplied_wall_hist = nan(N,1);
dQimplied_cap_hist  = nan(N,1);
dEimplied_hist      = nan(N,1);
residAb_hist        = nan(N,1);

% ---- comparison heat histories ----
dQnet_masked_hist   = nan(N,1);
dQnet_unmasked_hist = nan(N,1);
dqnet_Tnew_hist     = nan(N,1);  % inner evaluated with Tnew
dqnet_Tprev_hist    = nan(N,1);  % inner evaluated with Tprev

% store initial energy (t=0 state)
E_solid_hist(1) = solid_internal_energy(T, r_faces, dtheta, dz, rhoCp_ring) ...
    + endcap_internal_energy(Tcap0, TcapL, r_faces, cfg, BC);

snap = struct('t',{},'r_centers',{},'theta',{},'T',{});

% ---- Step 2: axial inner-wall insulation mask ----
useWallInsul = isfield(cfg,'wall_insul') && isfield(cfg.wall_insul,'enable') && cfg.wall_insul.enable;

if useWallInsul
    z_centers = ((1:Nz) - 0.5) * dz;
    m_in = ones(Nz,1);  % multiplier for inner wall gas BC

    if strcmpi(cfg.wall_insul.end,'z0')
        mask = z_centers <= cfg.wall_insul.length;
    elseif strcmpi(cfg.wall_insul.end,'zL')
        mask = z_centers >= (max(z_centers) - cfg.wall_insul.length);
    else
        error('cfg.wall_insul.end must be ''z0'' or ''zL''.');
    end

    m_in(mask) = cfg.wall_insul.factor;
else
    m_in = ones(Nz,1);
end

% Storage buffers
Tout3D_store = [];
t_store      = [];

%% =========================
% MAIN TIME LOOP
%% =========================
cache_rebuild_count = 0;
for step = 0:nsteps
    t = step * time.dt;
    t_step_total = tic;
    % ---- per-step timing locals (so "Other" is a true residual) ----
    store_time = 0; cache_time = 0; asm_time = 0; prec_time = 0; solve_time = 0;
    diagAb_time = 0; bc_time = 0; del_time = 0;

    % ---- store ----
    t_store_blk = tic;
    if mod(step, cfg.store_every) == 0
        T2 = mean(T,3);  % axial average for plotting/stats

        snap(end+1).t = t;
        snap(end).r_centers = r_centers;
        snap(end).theta = theta_centers;
        snap(end).T = T2;

        t_store(end+1,1) = t;
        Tout_theta_z = squeeze(T(end,:,:));        % [Ntheta x Nz]
        Tout3D_store(:,:,end+1) = Tout_theta_z;    % append along 3rd dim

        snap(end).Tbar_r = mean(T2, 2);
        snap(end).Tmax_r = max(T2, [], 2);
    end
    store_time = toc(t_store_blk);
    timing.store = timing.store + store_time;


    if step == nsteps
        break;
    end


    %% New Func
    [isBurning,T_gas,ndel,willDelete,willDel_hist,Lg,grain_present,n_grain] = compute_step_state(step, t, T, r_centers, dz, burn, geom, BC, cfg);

    % if mod(step, cfg.status_every) == 0
    %     fprintf('Lg = %.4f m | grain slices = %d / %d\n', Lg, n_grain, Nz);
    % end

    % ---- assemble + solve (cached base matrix; per-step BC updates) ----
    Tprev = T;
    t_start_step = tic;

    Tcap0_prev = Tcap0;
    TcapL_prev = TcapL;

    Nr     = numel(r_centers);
    Ntheta = size(Tprev,2);
    Nz     = size(Tprev,3);
    Nwall  = Nr*Ntheta*Nz;


    % Build cache once (or rebuild after any deletion that changes Nr / r_faces / r_centers)
    t_cache = tic;

    rebuildCache = ~exist('cache3D','var') || isempty(cache3D) || ...
        cache3D.Nr ~= numel(r_centers) || ...
        cache3D.Ntheta ~= size(Tprev,2) || ...
        cache3D.Nz ~= size(Tprev,3);
    
    if rebuildCache
        cache_rebuild_count = cache_rebuild_count + 1;
    
        cache3D = build_BE_polar3D_cache( ...
            r_faces, r_centers, dtheta, dz, ...
            k_ring, rhoCp_ring, time.dt, cfg, Lg, size(Tprev,2), size(Tprev,3));
    end
    if rebuildCache && mod(step, cfg.status_every) == 0
        fprintf('   cache rebuilt #%d (%.3f s)\n', cache_rebuild_count, cache_time);
    end
    if ~isfield(cache3D,'Mdiag') || isempty(cache3D.Mdiag)
        cache3D.Mdiag = build_Mdiag_wall_endcap(r_faces, r_centers, dtheta, dz, rhoCp_ring, time.dt, cfg);
    end

    cache_time = toc(t_cache);
    timing.cache = timing.cache + cache_time;


    t_asm = tic;

    [A, b] = assemble_from_cache_BE_polar3D( ...
        cache3D, Tprev, Tcap0_prev, TcapL_prev, BC, T_gas, cfg, Lg);
    asm_time = toc(t_asm);

    t_solve = tic;
    timing.assembly = timing.assembly + asm_time;
    % ---- endcap presence flags (must exist before pack_state_vec) ----
    cap_z0 = false;
    cap_zL = false;
    
    if isfield(cfg,'endcap') && isfield(cfg.endcap,'enable') && cfg.endcap.enable
        if ~isfield(cfg.endcap,'ends'); cfg.endcap.ends = 'both'; end
    
        cap_z0 = strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'z0');
        cap_zL = strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'zL');
    end


    %% =========================
    % Solve with rails: iterative with checks + retry + direct fallback
    %% =========================

    % --- Ensure cfg.solve exists ---
    if ~isfield(cfg,'solve'); cfg.solve = struct(); end

    % --- Diagnostic mode flag controls solver strictness ---
    diagMode = false;
    if isfield(cfg,'doEnergyDiagnostics') && cfg.doEnergyDiagnostics
        diagMode = true;
    end


    % Defaults (you can tune)
    if ~isfield(cfg.solve,'tol_diag');    cfg.solve.tol_diag = 1e-12; end
    if ~isfield(cfg.solve,'tol_prod');    cfg.solve.tol_prod = 1e-8;  end
    if ~isfield(cfg.solve,'maxit_diag');  cfg.solve.maxit_diag = 400; end
    if ~isfield(cfg.solve,'maxit_prod');  cfg.solve.maxit_prod = 200; end

    if ~isfield(cfg.solve,'sys_resid_max_diag'); cfg.solve.sys_resid_max_diag = 1e-11; end
    if ~isfield(cfg.solve,'sys_resid_max_prod'); cfg.solve.sys_resid_max_prod = 1e-8;  end

    tol   = diagMode * cfg.solve.tol_diag   + (~diagMode) * cfg.solve.tol_prod;
    maxit = diagMode * cfg.solve.maxit_diag + (~diagMode) * cfg.solve.maxit_prod;
    sys_resid_max = diagMode * cfg.solve.sys_resid_max_diag + (~diagMode) * cfg.solve.sys_resid_max_prod;

    % --- Warm start vector ---
    x0 = pack_state_vec(Tprev, cap_z0, cap_zL, Tcap0_prev, TcapL_prev);

    % --- Decide if we must rebuild preconditioner ---
    N = size(A,1);
    needPrec = true;
    if isfield(cache3D,'precN') && cache3D.precN == N
        needPrec = false;
    end
    if exist('didDeleteThisStep','var') && didDeleteThisStep
        needPrec = true;
    end
    if ~isfield(cache3D,'precRefreshEvery'); cache3D.precRefreshEvery = 50; end
    if ~isfield(cache3D,'stepCount'); cache3D.stepCount = 0; end
    cache3D.stepCount = cache3D.stepCount + 1;
    if mod(cache3D.stepCount, cache3D.precRefreshEvery) == 0
        needPrec = true;
    end

    % --- Build preconditioner if needed ---
    t_prec = tic;
    if needPrec
        if ~isfield(cfg.solve,'ichol_droptol');   cfg.solve.ichol_droptol = 1e-3; end
        if ~isfield(cfg.solve,'ichol_diagcomp');  cfg.solve.ichol_diagcomp = 1e-3; end
        if ~isfield(cfg.solve,'ilu_droptol');     cfg.solve.ilu_droptol = 1e-3; end

        cache3D.precOK = false;
        try
            cache3D.M = ichol(A, struct('type','ict', ...
                'droptol', cfg.solve.ichol_droptol, ...
                'diagcomp', cfg.solve.ichol_diagcomp));
            cache3D.Mt = cache3D.M';
            cache3D.precOK = true;
        catch
            setup = struct('type','ilutp','droptol',cfg.solve.ilu_droptol);
            [cache3D.L, cache3D.U] = ilu(A, setup);
            cache3D.precOK = false; % gmres path
        end
        cache3D.precN = N;
    end
    prec_time = toc(t_prec);
    timing.prec = timing.prec + prec_time;

    % --- Solve with rails (retry tightening then direct fallback) ---
    Tvec = [];
    solverName = '';
    flag = 0; relres = NaN; iters = [];

    for attempt = 1:3
        if isfield(cache3D,'precOK') && cache3D.precOK
            solverName = 'pcg';
            [Tvec, flag, relres, iters] = pcg(A, b, tol, maxit, cache3D.M, cache3D.Mt, x0);
        else
            solverName = 'gmres';
            restart = 75;
            [Tvec, flag, relres, iters] = gmres(A, b, restart, tol, maxit, cache3D.L, cache3D.U, x0);
        end

        % System residual check
        denom = max(norm(b), 1.0);
        sys_resid = norm(A*Tvec - b) / denom;

        % Accept if both iterative convergence and system residual are good
        if (flag == 0) && (sys_resid <= sys_resid_max)
            break;
        end

        % If not accepted, tighten and retry
        tol = tol * 0.1;
        maxit = maxit * 2;

        % On the second failure, force a preconditioner refresh next time
        if attempt == 2
            needPrec = true;
            cache3D.precN = -1; % guarantees rebuild on next step if you keep logic
        end
    end

    % Final safety fallback: direct solve if still not acceptable
    denom = max(norm(b), 1.0);
    sys_resid = norm(A*Tvec - b) / denom;
    if ~(flag == 0 && sys_resid <= sys_resid_max)
        if diagMode
            warning('Iterative solve not accepted (flag=%d, sys_resid=%.3g, relres=%.3g). Falling back to direct A\\b.', ...
                flag, sys_resid, relres);
        end
        solverName = 'direct';
        Tvec = A \ b;
        sys_resid = norm(A*Tvec - b) / denom;
    end

    solve_time = toc(t_solve);
    timing.solve = timing.solve + solve_time;

    % Unpack solution ONCE (no duplicate reshapes)
    [T, Tcap0, TcapL] = unpack_state_vec(Tvec, Nr, Ntheta, Nz, cap_z0, cap_zL);

    % Optional system residual check (cheap-ish; keep gated)
    if isfield(cfg,'diag') && isfield(cfg.diag,'consistency') && cfg.diag.consistency ...
            && mod(step, cfg.status_every) == 0
        sys_resid = norm(A*Tvec - b) / max(norm(b),1);
    end

    % ---- step time so far ----
    step_time = toc(t_start_step);

    %% =========================
    % Ab-implied energy and heat diagnostic (pre-deletion)
    %% =========================
    
    Tprev_vec = Tprev(:);
    if isfield(cfg,'endcap') && isfield(cfg.endcap,'enable') && cfg.endcap.enable
        if ~isfield(cfg.endcap,'ends'); cfg.endcap.ends = 'both'; end
        if strcmpi(cfg.endcap.ends,'both')
            Tprev_vec = [Tprev_vec; Tcap0_prev; TcapL_prev];
        elseif strcmpi(cfg.endcap.ends,'z0')
            Tprev_vec = [Tprev_vec; Tcap0_prev];
        elseif strcmpi(cfg.endcap.ends,'zL')
            Tprev_vec = [Tprev_vec; TcapL_prev];
        end
    end

    t_diagAb = tic;
    diagAb = energy_diag_from_Ab(A, b, Tvec, Tprev_vec, cache3D.Mdiag, time.dt, Nwall);
    diagAb_time = toc(t_diagAb);
    timing.diagAb = timing.diagAb + diagAb_time;


    dQimplied_hist(step+1)      = diagAb.dQ_implied_total;
    dQimplied_wall_hist(step+1) = diagAb.dQ_implied_wall;
    dQimplied_cap_hist(step+1)  = diagAb.dQ_implied_cap;
    dEimplied_hist(step+1)      = diagAb.dE_mass_total;
    residAb_hist(step+1)        = diagAb.resid_total;

    %% =========================
    % Boundary heat diagnostics  (kept aligned with your existing logic)
    %% =========================
    dQ_in_unmasked_Tnew  = 0.0;
    dQ_in_unmasked_Tprev = 0.0;

    dQ_in_unmasked  = 0.0;
    dQ_in_masked    = 0.0;
    dQ_out          = 0.0;   % positive leaving solid

    zmask_net = grain_present(:);   %#ok<NASGU>

    qface_inner_net      = zeros(Ntheta, Nz);  %#ok<NASGU>
    qface_inner_unmasked = zeros(Ntheta, Nz);  %#ok<NASGU>

    n_inGrain_count = 0;

    t_bc = tic;
    for kk = 1:Nz
        inGrain = grain_present(kk);
        n_inGrain_count = n_inGrain_count + double(inGrain);

        if ~isempty(m_in)
            m = m_in(kk);
        else
            m = 1.0;
        end

        for j = 1:Ntheta
            % INNER boundary (i=1)
            Tw0 = Tprev(1,j,kk);
            [h_eff, q_const] = rad_lin_coeff(BC.eps_in, BC.sigma, Tw0, T_gas);
            h_base = (BC.h_in + h_eff);

            Aface_in = r_faces(1) * dtheta * dz;

            dQ_face_in_Tnew  = time.dt * ( h_base*Aface_in*(T_gas - T(1,j,kk))      + q_const*Aface_in );
            dQ_face_in_Tprev = time.dt * ( h_base*Aface_in*(T_gas - Tprev(1,j,kk)) + q_const*Aface_in );

            dQ_in_unmasked_Tnew  = dQ_in_unmasked_Tnew  + dQ_face_in_Tnew;
            dQ_in_unmasked_Tprev = dQ_in_unmasked_Tprev + dQ_face_in_Tprev;

            dQ_in_unmasked = dQ_in_unmasked + dQ_face_in_Tnew;
            qface_inner_unmasked(j,kk) = dQ_face_in_Tnew;

            if inGrain
                h_total = m * h_base;
                dQ_face_in_masked = time.dt * ( h_total*Aface_in*(T_gas - T(1,j,kk)) + q_const*Aface_in );
                dQ_in_masked = dQ_in_masked + dQ_face_in_masked;
                qface_inner_net(j,kk) = dQ_face_in_masked;
            else
                qface_inner_net(j,kk) = 0.0;
            end

            % OUTER boundary (i=Nr)
            Tw1 = Tprev(Nr,j,kk);
            [h_eff_o, q_const_o] = rad_lin_coeff(BC.eps_out, BC.sigma, Tw1, BC.T_amb);
            h_total_o = BC.h_out + h_eff_o;

            Aface_out = r_faces(end) * dtheta * dz;
            dQ_into_solid_outer = time.dt * ( h_total_o*Aface_out*(BC.T_amb - T(Nr,j,kk)) + q_const_o*Aface_out );

            % dQ_out positive leaving solid
            dQ_out = dQ_out - dQ_into_solid_outer;
        end
    end
    bc_time = toc(t_bc);
    timing.bcHeat = timing.bcHeat + bc_time;


    dQnet_masked   = dQ_in_masked   - dQ_out;
    dQnet_unmasked = dQ_in_unmasked - dQ_out;

    dqnet_Tnew  = dQ_in_unmasked_Tnew  - dQ_out;
    dqnet_Tprev = dQ_in_unmasked_Tprev - dQ_out;

    dQnet_masked_hist(step+1)   = dQnet_masked;
    dQnet_unmasked_hist(step+1) = dQnet_unmasked;
    dqnet_Tnew_hist(step+1)     = dqnet_Tnew;
    dqnet_Tprev_hist(step+1)    = dqnet_Tprev;

    % Use Ab-implied heat for closure while diagnosing, as you chose
    dQnet = diagAb.dQ_implied_total;

    %% =========================
    % Energy balance for this step (pre-deletion state)
    %% =========================
    E_prev = E_solid_hist(step+1);

    E_wall_new = solid_internal_energy(T, r_faces, dtheta, dz, rhoCp_ring);
    E_cap_new  = endcap_internal_energy(Tcap0, TcapL, r_faces, cfg, BC);
    E_new      = E_wall_new + E_cap_new;

    % Store pre-deletion new energy first
    E_solid_hist(step+2) = E_new;

    %% =========================
    % Energy removed if we delete ndel rings after solve
    %% =========================
    E_removed_wall = 0.0;

    t_del = tic;
    if willDelete
        for ii = 1:ndel
            r_imh = r_faces(ii);
            r_iph = r_faces(ii+1);
            Vcell = 0.5*(r_iph^2 - r_imh^2) * dtheta * dz;
            T_ring = squeeze(T(ii,:,:));
            E_removed_wall = E_removed_wall + rhoCp_ring(ii) * Vcell * sum(T_ring(:));
        end
    end
    del_time = toc(t_del);
    timing.deletion = timing.deletion + del_time;

    Erem_hist(step+1) = E_removed_wall;

    E_removed_cap = 0.0;
    if willDelete && isfield(cfg,'endcap') && isfield(cfg.endcap,'enable') && cfg.endcap.enable
        if ~isfield(cfg.endcap,'ends'); cfg.endcap.ends = 'both'; end

        rhoCp_cap = cfg.endcap.rho * cfg.endcap.cp;

        r_outer = r_faces(end);
        r_inner_before = r_faces(1);
        r_inner_after  = r_faces(ndel+1);

        A_before = pi * (r_outer^2 - r_inner_before^2);
        A_after  = pi * (r_outer^2 - r_inner_after^2);
        dV = (A_before - A_after) * cfg.endcap.thickness;

        if dV > 0
            if strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'z0')
                if ~isempty(Tcap0)
                    E_removed_cap = E_removed_cap + rhoCp_cap * dV * Tcap0;
                end
            end
            if strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'zL')
                if ~isempty(TcapL)
                    E_removed_cap = E_removed_cap + rhoCp_cap * dV * TcapL;
                end
            end
        end
    end

    E_removed_total = E_removed_wall + E_removed_cap;

    %% =========================
    % Closure defined on the solved domain (pre-deletion)
    %% =========================
    dE_corr = (E_new - E_prev);
    closure = dE_corr - dQnet;

    Qscale = max(abs(dQ_in_masked) + abs(dQ_out), 1.0);
    relerr = closure / Qscale;

    dE_hist(step+1)     = dE_corr;
    dQnet_hist(step+1)  = dQnet;
    relerr_hist(step+1) = relerr;

    closure_pre_hist(step+1) = closure;

    %% =========================
    % Optional diagnostic: conduction across the DELETION INTERFACE from A
    %% =========================
    Qcond12_fromA = NaN;

    if willDelete
        i_del  = ndel;
        i_keep = ndel + 1;

        if i_del >= 1 && i_keep <= Nr
            [Jgrid, Kgrid] = ndgrid(1:Ntheta, 1:Nz);
            Jlist = Jgrid(:);
            Klist = Kgrid(:);

            p_del  = sub2ind([Nr, Ntheta, Nz], i_del *ones(size(Jlist)),  Jlist, Klist);
            p_keep = sub2ind([Nr, Ntheta, Nz], i_keep*ones(size(Jlist)),  Jlist, Klist);

            lin = sub2ind(size(A), p_del, p_keep);
            Gface = full(-A(lin));

            dT = reshape(T(i_del,:,:) - T(i_keep,:,:), [], 1);

            Qcond12_fromA = time.dt * sum(Gface .* dT);
        end
    end

    Qcond12A_hist(step+1) = Qcond12_fromA;

    %% ---- status output ----
    if mod(step, cfg.status_every) == 0 || step == 1
        Tmax_now = max(T(:));
        Tmin_now = min(T(:));
        Tout_max_now = max(T(end,:,:),[],'all');

        if isBurning
            r_burned = geom.r_inner0 + burn.rburn * t;
            burn_msg = sprintf('burning, r_front=%.4f m', r_burned);
        else
            burn_msg = 'post-burn';
        end

        fprintf(['step %5d / %5d | t = %6.2f s | %s | ' ...
            'Nr = %4d | Nθ = %3d | cells = %7d | ' ...
            'T[min,max] = [%6.1f,%6.1f] K | ' ...
            'T_out,max = %6.1f K | ' ...
            'asm = %.2fs | solve = %.2fs | total = %.2fs\n'], ...
            step, nsteps, t, burn_msg, ...
            Nr, Ntheta, Nr*Ntheta*Nz, ...
            Tmin_now, Tmax_now, Tout_max_now, ...
            asm_time, solve_time, step_time);

        if step_time > cfg.time_warn
            fprintf('   step time exceeded %.1f s\n', cfg.time_warn);
        end
    end

    %% =========================
    % Deletion happens after diagnostics, then set carried energy correctly
    %% =========================
    if willDelete
        Nr_before = numel(r_centers);

        % perform deletion on the solved temperature field
        T(1:ndel,:,:)             = [];
        r_centers(1:ndel)         = [];
        r_faces(1:ndel)           = [];
        layer_idx_of_ring(1:ndel) = [];

        Nr_after = numel(r_centers);
        if exist('cache3D','var') && isstruct(cache3D)
            cache3D.Mdiag = [];
        end

        if Nr_after ~= (Nr_before - ndel)
            error('Nr observed mismatch: before=%d after=%d ndel=%d', Nr_before, Nr_after, ndel);
        end

        % Update per-ring properties after deletion
        [k_ring, rhoCp_ring] = ring_props(layers, layer_idx_of_ring);

        % Energies after deletion
        E_afterDel_wall = solid_internal_energy(T, r_faces, dtheta, dz, rhoCp_ring);
        E_afterDel_cap  = endcap_internal_energy(Tcap0, TcapL, r_faces, cfg, BC);
        E_afterDel      = E_afterDel_wall + E_afterDel_cap;

        % Overwrite carried energy with post-deletion energy (this is crucial)
        E_solid_hist(step+2) = E_afterDel;

        Qdel = E_removed_total;
        Qdel_hist(step+1) = Qdel;

        R_post = (E_afterDel - E_prev) - dQnet + Qdel;
        closure_post_hist(step+1) = R_post;
        closure_hist(step+1)      = R_post;

    else
        Qdel_hist(step+1) = 0.0;
        closure_post_hist(step+1) = closure_pre_hist(step+1);
        closure_hist(step+1)      = closure_pre_hist(step+1);
    end
    
    % Instead, do the robust way: use the timing buckets you already accumulate.
    % The per-step total is step_total; you subtract the per-step contributions you just added.
    % Easiest drop-in: compute other as a residual from step_total and the known block durations you measured this step.
    
    % ---- drop-in: track per-step block durations explicitly ----
    % Right after each toc(...) assign to a local var and add to timing.
    % Then:
    step_total = toc(t_step_total);
    known_step = store_time + cache_time + asm_time + prec_time + solve_time + diagAb_time + bc_time + del_time;
    timing.other = timing.other + max(step_total - known_step, 0);


end

total = sum(struct2array(timing));

fprintf('\n===== TIMING SUMMARY =====\n');
fprintf('Total runtime        : %.2f s\n', total);
fprintf('Cache rebuild        : %.2f s (%.1f%%)\n', timing.cache,    100*timing.cache/total);
fprintf('Assembly             : %.2f s (%.1f%%)\n', timing.assembly, 100*timing.assembly/total);
fprintf('Preconditioner build : %.2f s (%.1f%%)\n', timing.prec,     100*timing.prec/total);
fprintf('Linear solve         : %.2f s (%.1f%%)\n', timing.solve,    100*timing.solve/total);
fprintf('Ab diagnostics       : %.2f s (%.1f%%)\n', timing.diagAb,   100*timing.diagAb/total);
fprintf('BC heat diagnostics  : %.2f s (%.1f%%)\n', timing.bcHeat,   100*timing.bcHeat/total);
fprintf('Deletion bookkeeping : %.2f s (%.1f%%)\n', timing.deletion, 100*timing.deletion/total);
fprintf('Other                : %.2f s (%.1f%%)\n', timing.other,    100*timing.other/total);
fprintf('==========================\n');

%% Final axial profiles
z_centers = ((1:cfg.Nz) - 0.5) * dz;
Tin_z = squeeze(mean(T(1,:,:), 2));
Tout_z = squeeze(mean(T(end,:,:), 2));

%% =========================================================
% FINAL OUTER-WALL TEMPERATURE SUMMARY
% =========================================================

% Safety checks
assert(exist('Tout3D_store','var') == 1, 'Tout3D_store not found');
assert(exist('t_store','var') == 1, 't_store not found');
assert(exist('z_centers','var') == 1, 'z_centers not found');

Nt = size(Tout3D_store, 3);

%% -------- 1) Hottest time at a specific z-slice --------
% Choose the z-slice you already plot or care about
% Examples:
%   z_query = Lg;              % grain end
%   z_query = 0.0;             % bulkhead
%   z_query = max(z_centers);  % nozzle end

z_query = Lg;   % <<< CHANGE IF DESIRED

[~, iz] = min(abs(z_centers - z_query));

% Collapse theta, keep time
Tz_time = squeeze(max(Tout3D_store(:, iz, :), [], 1));  % [Nt x 1]

[Tz_max, it_zmax] = max(Tz_time);
t_zmax = t_store(it_zmax);

fprintf('\n===== OUTER WALL @ z = %.4f m =====\n', z_centers(iz));
fprintf('Max temperature : %.2f K\n', Tz_max);
fprintf('Time of max     : %.2f s\n', t_zmax);

%% -------- 2) Global hottest outer-wall point (entire run) --------
[Tglobal_max, idx] = max(Tout3D_store(:));

[jmax, kmax, tmax] = ind2sub(size(Tout3D_store), idx);

fprintf('\n===== GLOBAL OUTER WALL MAX =====\n');
fprintf('Max temperature : %.2f K\n', Tglobal_max);
fprintf('Time            : %.2f s\n', t_store(tmax));
fprintf('Axial location  : z = %.4f m\n', z_centers(kmax));
fprintf('Theta index     : %d / %d\n', jmax, size(Tout3D_store,1));
fprintf('=================================\n\n');


figure;
plot(z_centers, Tin_z, 'LineWidth', 2); hold on;
plot(z_centers, Tout_z, 'LineWidth', 2);
grid on;
xlabel('z [m]');
ylabel('Temperature [K]');
legend('Inner wall (mean over \theta)','Outer wall (mean over \theta)','Location','best');
title('Axial temperature profiles at final time');

plot_results(snap, layers, geom, burn, BC);

%% ==================== 3D tube animation + GIF export (outer surface) ====================
if exist('Tout3D_store','var') && ~isempty(Tout3D_store)
    r_outer = r_faces(end);
    z_centers = ((1:cfg.Nz) - 0.5) * (cfg.L / cfg.Nz);

    [TH, ZZ] = meshgrid(theta_centers, z_centers);
    XX = r_outer * cos(TH);
    YY = r_outer * sin(TH);

    C0 = Tout3D_store(:,:,1).';   % [Nz x Ntheta]

    fig = figure('Name','3D Outer Surface Temperature', 'Color','w');
    h = surf(XX, YY, ZZ, C0, 'EdgeColor','none');
    axis equal tight;
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    cb = colorbar; ylabel(cb, 'Temperature [K]');
    colormap(turbo);
    view(35, 20);
    camlight headlight; lighting gouraud;

    Tmin = min(Tout3D_store, [], 'all');
    Tmax = max(Tout3D_store, [], 'all');
    caxis([Tmin Tmax]);

    NtT = size(Tout3D_store, 3);
    Ntt = numel(t_store);
    Nt  = min(NtT, Ntt);

    gifname = 'outer_wall_temperature.gif';
    gif_dt  = 0.05;

    for kk = 1:Nt
        Ck = Tout3D_store(:,:,kk).';
        set(h, 'CData', Ck);
        title(sprintf('Outer wall T, t = %.2f s', t_store(kk)));
        drawnow;

        frame = getframe(fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        if kk == 1
            imwrite(imind, cm, gifname, 'gif', 'Loopcount', inf, 'DelayTime', gif_dt);
        else
            imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime', gif_dt);
        end
    end

    fprintf('GIF written to %s (%d frames)\n', gifname, Nt);
else
    warning('Tout3D_store not found or empty. Add the store block inside the sim loop first.');
end

%% Energy diagnostic plots
if isfield(cfg,'doEnergyDiagnostics') && cfg.doEnergyDiagnostics
    cfg2 = cfg;

    cfg2.diagLabel = 'POST (after deletion correction)';
    plot_energy_diagnostics(time, closure_post_hist, relerr_hist, willDel_hist, Qcond12A_hist, cfg2);
end

end

%% ==================== Helper functions ====================

function [dtheta, theta_centers] = uniform_theta_grid(Ntheta)
dtheta = 2*pi / Ntheta;
theta_centers = ( (0:Ntheta-1) + 0.5 ) * dtheta;
end

function [r_faces, r_centers, layer_idx_of_ring] = build_radial_faces_with_interfaces(geom, layers, dr)
% Build radial faces from r_inner0 to r_outer with interfaces on faces
r0 = geom.r_inner0;

total_thick = 0.0;
for i = 1:numel(layers)
    total_thick = total_thick + layers(i).thickness;
end
%r_outer = r0 + total_thick;

% Build interface radii (faces)
interfaces = zeros(numel(layers)+1,1);
interfaces(1) = r0;
acc = r0;
for i = 1:numel(layers)
    acc = acc + layers(i).thickness;
    interfaces(i+1) = acc;
end

% Build faces ensuring each interface is a face location
r_faces = r0;
%layer_idx_of_ring = [];

for li = 1:numel(layers)
    rL = interfaces(li);
    rR = interfaces(li+1);

    % Ensure current face is exactly rL
    if abs(r_faces(end) - rL) > 1e-12
        r_faces(end+1,1) = rL;
    end

    % March by dr until reaching rR
    while true
        r_next = r_faces(end) + dr;
        if r_next >= rR - 1e-12
            if abs(r_faces(end) - rR) > 1e-12
                r_faces(end+1,1) = rR;
            end
            break;
        else
            r_faces(end+1,1) = r_next;
        end
    end
end

% Build centers and assign each ring to a layer by midpoint radius
Nr = numel(r_faces) - 1;
r_centers = 0.5*(r_faces(1:Nr) + r_faces(2:Nr+1));

layer_idx_of_ring = zeros(Nr,1);
for i = 1:Nr
    rmid = r_centers(i);
    % find which interface interval contains rmid
    for li = 1:numel(layers)
        if (rmid >= interfaces(li) - 1e-12) && (rmid <= interfaces(li+1) + 1e-12)
            layer_idx_of_ring(i) = li;
            break;
        end
    end
end
end

function [k_ring, rhoCp_ring] = ring_props(layers, layer_idx_of_ring)
Nr = numel(layer_idx_of_ring);
k_ring = zeros(Nr,1);
rhoCp_ring = zeros(Nr,1);
for i = 1:Nr
    li = layer_idx_of_ring(i);
    mat = layers(li).material;
    k_ring(i) = mat.k;
    rhoCp_ring(i) = mat.rho * mat.cp;
end
end

function [A, b] = assemble_BE_polar_3D(Tprev, Tcap0_prev, TcapL_prev, r_faces, r_centers, dtheta, dz, ...
    k_ring, rhoCp_ring, dt, BC, ~, T_gas, cfg, Lg)
%ASSEMBLE_BE_POLAR_3D Backward Euler FV for 3D polar conduction (r,theta,z) with optional lumped endcaps.
% Unknown ordering for wall cells: p = sub2ind([Nr, Ntheta, Nz], i, j, k)
% Optional endcap nodes appended after wall cells: [Tcap0; TcapL] depending on cfg.endcap.ends.
%
% Step 1 integration: optionally disable hot-gas convection/radiation on an endcap inner face
% via cfg.endcap.gas_exposed_z0 and cfg.endcap.gas_exposed_zL (true/false).

Nr     = numel(r_centers);
Ntheta = size(Tprev,2);
Nz     = size(Tprev,3);
Nwall  = Nr * Ntheta * Nz;

useEndcap = isfield(cfg,'endcap') && isfield(cfg.endcap,'enable') && cfg.endcap.enable;

% ---- Step 4: compute throat-area proxy scaling for cavity convection ----
useNozProxy = isfield(cfg,'noz') && isfield(cfg.noz,'enableProxy') && cfg.noz.enableProxy;

% ---- Step 4: throat-area proxy scaling needs port area ----
r_inner = r_faces(1);           % current inner radius (after deletion)
Aport   = pi * r_inner^2;       % current port cross-sectional area


r_t = [];
At  = [];
if useNozProxy
    r_t = cfg.noz.r_t;
    At  = pi * r_t^2;                 % throat area
end

% For cavity-side "reference area", use current port cross-section area
Aport = pi * r_inner^2;               % uses current inner radius (after deletion)


% ---- Axial grain-present mask (grain shrinks from z=L backward) ----
z_centers = ((1:Nz) - 0.5) * dz;

% ---- Step 2: inner-wall insulation multiplier along z (default 1 everywhere) ----
m_in = ones(Nz,1);

useWallInsul = isfield(cfg,'wall_insul') && isfield(cfg.wall_insul,'enable') && cfg.wall_insul.enable;
if useWallInsul
    if ~isfield(cfg.wall_insul,'end');    cfg.wall_insul.end = 'z0'; end
    if ~isfield(cfg.wall_insul,'length'); cfg.wall_insul.length = 0.0; end
    if ~isfield(cfg.wall_insul,'factor'); cfg.wall_insul.factor = 0.0; end

    if strcmpi(cfg.wall_insul.end,'z0')
        mask = (z_centers <= cfg.wall_insul.length);
    elseif strcmpi(cfg.wall_insul.end,'zL')
        mask = (z_centers >= (max(z_centers) - cfg.wall_insul.length));
    else
        error('cfg.wall_insul.end must be ''z0'' or ''zL''.');
    end

    m_in(mask) = cfg.wall_insul.factor;
end

% Grain occupies 0 <= z <= Lg. Aft cavity is z > Lg.
grain_present = (z_centers <= Lg + 1e-12);

cap_z0 = false;
cap_zL = false;
if useEndcap
    if ~isfield(cfg.endcap,'ends'); cfg.endcap.ends = 'both'; end
    cap_z0 = strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'z0');
    cap_zL = strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'zL');
end
Ncap = double(cap_z0) + double(cap_zL);
N    = Nwall + Ncap;

% Per-end gas exposure flags (defaults true)
gas_exposed_z0 = true;
gas_exposed_zL = true;
if useEndcap && isfield(cfg.endcap,'gas_exposed_z0'); gas_exposed_z0 = cfg.endcap.gas_exposed_z0; end
if useEndcap && isfield(cfg.endcap,'gas_exposed_zL'); gas_exposed_zL = cfg.endcap.gas_exposed_zL; end

% Preallocate sparse triplets
nz_est = N * 14;
I = zeros(nz_est,1);
J = zeros(nz_est,1);
V = zeros(nz_est,1);
b = zeros(N,1);
idx = 0;

% Periodic theta indexing
jp = [2:Ntheta 1];
jm = [Ntheta 1:Ntheta-1];


% Endcap node indices (if present)
p_cap0 = [];
p_capL = [];
if cap_z0
    p_cap0 = Nwall + 1;
end
if cap_zL
    p_capL = Nwall + 1 + double(cap_z0);
end

% Annulus area for endcap convection/radiation (matches wall domain)
r_inner   = r_faces(1);
r_outer   = r_faces(end);
A_annulus = pi * (r_outer^2 - r_inner^2);  % [m^2]

% Convenience: pull endcap parameters once (no repeated isfield inside loops)
if useEndcap
    if ~isfield(cfg.endcap,'thickness'); cfg.endcap.thickness = dz; end
    if ~isfield(cfg.endcap,'k');         cfg.endcap.k         = 1.0; end
    if ~isfield(cfg.endcap,'rho');       cfg.endcap.rho       = 1000; end
    if ~isfield(cfg.endcap,'cp');        cfg.endcap.cp        = 1000; end
    if ~isfield(cfg.endcap,'h_in');      cfg.endcap.h_in      = 0; end
    if ~isfield(cfg.endcap,'eps_in');    cfg.endcap.eps_in    = 0; end
    if ~isfield(cfg.endcap,'h_out');     cfg.endcap.h_out     = 0; end
    if ~isfield(cfg.endcap,'eps_out');   cfg.endcap.eps_out   = 0; end
end

% Assemble endcap node equations (lumped) if enabled
if useEndcap
    rhoCp_cap = cfg.endcap.rho * cfg.endcap.cp;
    Vcap      = A_annulus * cfg.endcap.thickness;  % [m^3]
    aPcap     = rhoCp_cap * Vcap / dt;             % [W/K] in BE form

    if cap_z0
        pcap     = p_cap0;
        Tcap_prev = Tcap0_prev;

        % transient
        idx=idx+1; I(idx)=pcap; J(idx)=pcap; V(idx)=V(idx) + aPcap;
        b(pcap) = b(pcap) + aPcap * Tcap_prev;

        % inside (to gas): only if gas reaches this end
        if gas_exposed_z0
            Tw0 = Tcap_prev;
            [h_rad_in, q_const_in] = rad_lin_coeff(cfg.endcap.eps_in, BC.sigma, Tw0, T_gas);
            htot_in = cfg.endcap.h_in + h_rad_in;
            if htot_in > 0
                idx=idx+1; I(idx)=pcap; J(idx)=pcap; V(idx)=V(idx) + htot_in*A_annulus;
                b(pcap) = b(pcap) + htot_in*A_annulus*T_gas + q_const_in*A_annulus;
            end
        end

        % outside (to ambient)
        Tw1 = Tcap_prev;
        [h_rad_out, q_const_out] = rad_lin_coeff(cfg.endcap.eps_out, BC.sigma, Tw1, BC.T_amb);
        htot_out = cfg.endcap.h_out + h_rad_out;
        if htot_out > 0
            idx=idx+1; I(idx)=pcap; J(idx)=pcap; V(idx)=V(idx) + htot_out*A_annulus;
            b(pcap) = b(pcap) + htot_out*A_annulus*BC.T_amb + q_const_out*A_annulus;
        end
    end

    if cap_zL
        pcap     = p_capL;
        Tcap_prev = TcapL_prev;

        % transient
        idx=idx+1; I(idx)=pcap; J(idx)=pcap; V(idx)=V(idx) + aPcap;
        b(pcap) = b(pcap) + aPcap * Tcap_prev;

        % inside (to gas): only if gas reaches this end
        if gas_exposed_zL
            Tw0 = Tcap_prev;
            [h_rad_in, q_const_in] = rad_lin_coeff(cfg.endcap.eps_in, BC.sigma, Tw0, T_gas);
            htot_in = cfg.endcap.h_in + h_rad_in;
            if htot_in > 0
                idx=idx+1; I(idx)=pcap; J(idx)=pcap; V(idx)=V(idx) + htot_in*A_annulus;
                b(pcap) = b(pcap) + htot_in*A_annulus*T_gas + q_const_in*A_annulus;
            end
        end

        % outside (to ambient)
        Tw1 = Tcap_prev;
        [h_rad_out, q_const_out] = rad_lin_coeff(cfg.endcap.eps_out, BC.sigma, Tw1, BC.T_amb);
        htot_out = cfg.endcap.h_out + h_rad_out;
        if htot_out > 0
            idx=idx+1; I(idx)=pcap; J(idx)=pcap; V(idx)=V(idx) + htot_out*A_annulus;
            b(pcap) = b(pcap) + htot_out*A_annulus*BC.T_amb + q_const_out*A_annulus;
        end
    end
end

% Main wall assembly
for kk = 1:Nz
    for j = 1:Ntheta
        for i = 1:Nr
            p = sub2ind([Nr, Ntheta, Nz], i, j, kk);

            r_imh = r_faces(i);
            r_iph = r_faces(i+1);
            ri    = r_centers(i);

            % Control volume geometry (3D)
            Vcv = 0.5*(r_iph^2 - r_imh^2) * dtheta * dz;

            aP = rhoCp_ring(i) * Vcv / dt;
            bP = aP * Tprev(i,j,kk);

            % Radial conduction neighbors
            % ---------- Radial conduction / inner+outer BC ----------
            if i > 1
                % conduction to i-1
                r_face = r_imh;
                Aface  = r_face * dtheta * dz;
                dW     = ri - r_centers(i-1);
                k_face = harmonic_mean(k_ring(i-1), k_ring(i));
                GW     = k_face * Aface / dW;

                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + GW;
                q = p_all(i-1,j,kk);
                idx=idx+1; I(idx)=p; J(idx)=q; V(idx)=V(idx) - GW;

            else
                % i == 1: inner boundary
                r_face = r_imh;
                Aface  = r_face * dtheta * dz;

                if grain_present(kk)
                    % ---- Port region (grain present): use main hot-gas BC + Step 2 multiplier ----
                    m = m_in(kk);

                    Tw0 = Tprev(i,j,kk);
                    [h_eff, q_const] = rad_lin_coeff(BC.eps_in, BC.sigma, Tw0, T_gas);

                    h_total = m * (BC.h_in + h_eff);

                    if h_total > 0
                        idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + h_total*Aface;
                        bP = bP + h_total*Aface*T_gas + m*q_const*Aface;
                    end

                else
                    % ---- Aft cavity region (no grain): apply reduced gas BC ----
                    % Choose cavity gas temperature
                    if isfield(BC,'T_gas_cav_mode') && strcmpi(BC.T_gas_cav_mode,'ambient')
                        Tenv = BC.T_amb;
                    else
                        Tenv = T_gas;  % default: same as chamber gas for Step 3
                    end

                    % Get cavity parameters (with safe defaults)
                    % Get cavity emissivity
                    eps_cav = BC.eps_in_cav;

                    % Base cavity h
                    h_cav = cfg.noz.h_cav_base;

                    % Apply throat-area proxy scaling (optional)
                    if useNozProxy
                        % Ratio that increases when throat is "large" relative to port cross-section
                        % The exponent alpha lets you tune sensitivity without changing sign/behavior.
                        ratio = cfg.noz.CdA * (At / max(Aport, 1e-12));

                        h_cav = cfg.noz.h_cav_base * (ratio ^ cfg.noz.alpha);

                        % Clamp for numerical/physical sanity
                        h_cav = min(max(h_cav, cfg.noz.h_cav_min), cfg.noz.h_cav_max);
                    end


                    Tw0 = Tprev(i,j,kk);
                    [h_eff, q_const] = rad_lin_coeff(eps_cav, BC.sigma, Tw0, Tenv);

                    h_total = h_cav + h_eff;

                    if h_total > 0
                        idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + h_total*Aface;
                        bP = bP + h_total*Aface*Tenv + q_const*Aface;
                    end
                end
            end


            if i < Nr
                % conduction to i+1
                r_face = r_iph;
                Aface  = r_face * dtheta * dz;
                dE     = r_centers(i+1) - ri;
                k_face = harmonic_mean(k_ring(i), k_ring(i+1));
                GE     = k_face * Aface / dE;

                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + GE;
                q = p_all(i+1, j, kk);
                idx=idx+1; I(idx)=p; J(idx)=q; V(idx)=V(idx) - GE;

            else
                % i == Nr: outer boundary
                Tw0 = Tprev(i,j,kk);
                [h_eff, q_const] = rad_lin_coeff(BC.eps_out, BC.sigma, Tw0, BC.T_amb);
                h_total = BC.h_out + h_eff;

                r_face = r_iph;
                Aface  = r_face * dtheta * dz;

                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + h_total*Aface;
                bP = bP + h_total*Aface*BC.T_amb + q_const*Aface;
            end


            % Theta conduction neighbors (periodic)
            dr_i   = (r_iph - r_imh);
            Atheta = dr_i * dz;
            ds     = ri * dtheta;
            Gth    = k_ring(i) * Atheta / ds;

            jE = jp(j); jW = jm(j);
            qE = p_all(i, jE, kk);
            qW = p_all(i, jW, kk);

            idx=idx+1; I(idx)=p; J(idx)=p;  V(idx)=V(idx) + 2*Gth;
            idx=idx+1; I(idx)=p; J(idx)=qE; V(idx)=V(idx) - Gth;
            idx=idx+1; I(idx)=p; J(idx)=qW; V(idx)=V(idx) - Gth;

            % Axial conduction neighbors
            Az = 0.5*(r_iph^2 - r_imh^2) * dtheta;
            Gz = k_ring(i) * Az / dz;

            if kk > 1
                q = p_all(i,j,kk-1);
                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + Gz;
                idx=idx+1; I(idx)=p; J(idx)=q; V(idx)=V(idx) - Gz;
            else
                % z=0 end
                if cap_z0
                    Gcap = cfg.endcap.k * Az / cfg.endcap.thickness;

                    idx=idx+1; I(idx)=p;      J(idx)=p;      V(idx)=V(idx) + Gcap;
                    idx=idx+1; I(idx)=p;      J(idx)=p_cap0; V(idx)=V(idx) - Gcap;

                    idx=idx+1; I(idx)=p_cap0; J(idx)=p_cap0; V(idx)=V(idx) + Gcap;
                    idx=idx+1; I(idx)=p_cap0; J(idx)=p;      V(idx)=V(idx) - Gcap;
                else
                    if isfield(cfg,'endBC') && strcmpi(cfg.endBC,'robin')
                        h = BC.h_z0;
                        idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + h*Az;
                        bP = bP + h*Az*BC.T_amb;
                    end
                end
            end

            if kk < Nz
                q = p_all(i, j, kk+1);
                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + Gz;
                idx=idx+1; I(idx)=p; J(idx)=q; V(idx)=V(idx) - Gz;
            else
                % z=L end
                if cap_zL
                    Gcap = cfg.endcap.k * Az / cfg.endcap.thickness;

                    idx=idx+1; I(idx)=p;      J(idx)=p;      V(idx)=V(idx) + Gcap;
                    idx=idx+1; I(idx)=p;      J(idx)=p_capL; V(idx)=V(idx) - Gcap;

                    idx=idx+1; I(idx)=p_capL; J(idx)=p_capL; V(idx)=V(idx) + Gcap;
                    idx=idx+1; I(idx)=p_capL; J(idx)=p;      V(idx)=V(idx) - Gcap;
                else
                    if isfield(cfg,'endBC') && strcmpi(cfg.endBC,'robin')
                        h = BC.h_zL;
                        idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + h*Az;
                        bP = bP + h*Az*BC.T_amb;
                    end
                end
            end

            % Transient term
            idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + aP;
            b(p) = bP;
        end
    end
end

% Finalize sparse
I = I(1:idx);
J = J(1:idx);
V = V(1:idx);
A = sparse(I,J,V,N,N);

end

% function [Tavg, Tmax] = outer_wall_stats(T)
% Tout = T(end,:);
% Tavg = mean(Tout);
% Tmax = max(Tout);
% end

function plot_results(snap, layers, geom, burn, BC)
%PLOT_RESULTS Reproduce original 1D-style plots + 2D cross-section plot
% using theta-averaged profiles from the 2D solution.

if isempty(snap)
    error('plot_results: snap is empty.');
end

Nt = numel(snap);

% Time vector
t_vec = zeros(Nt,1);
for k = 1:Nt
    t_vec(k) = snap(k).t;
end

% ---- Outer surface temperature vs time (avg and max around theta) ----
Tout_avg = zeros(Nt,1);
Tout_max = zeros(Nt,1);
for k = 1:Nt
    Tout = snap(k).T(end,:);     % outer ring
    Tout_avg(k) = mean(Tout);
    Tout_max(k) = max(Tout);
end

figure;
plot(t_vec, Tout_avg, 'LineWidth', 2); hold on;
plot(t_vec, Tout_max, 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('Outer wall temperature [K]');
title('Outer surface temperature (avg and max around \theta)');
legend('avg', 'max', 'Location', 'best');


% ---- Build original-style radius-time matrices on a fixed reference radius ----
% Use the initial snapshot's radial centers as the fixed plotting grid
r_ref = snap(1).r_centers(:);
Nr_ref = numel(r_ref);
r_rel = r_ref - geom.r_inner0;   % match original x-axis definition

TbarMat = nan(Nt, Nr_ref);
TmaxMat = nan(Nt, Nr_ref);

for k = 1:Nt
    r_k = snap(k).r_centers(:);

    % theta-averaged and max profiles should exist
    if ~isfield(snap(k), 'Tbar_r')
        error('plot_results: snap(%d) missing Tbar_r. Store snap(end).Tbar_r = mean(T,2);', k);
    end
    Tbar_k = snap(k).Tbar_r(:);

    if isfield(snap(k), 'Tmax_r')
        Tmax_k = snap(k).Tmax_r(:);
    else
        % If not stored, compute now (slightly more expensive but OK for plotting)
        Tmax_k = max(snap(k).T, [], 2);
    end

    % Domain shrinks as burn removes inner rings. Treat missing inner region as "burned"
    r_min_k = min(r_k);

    % Fill burned region with gas temperature during burn, ambient after burn (to mimic your original viz)
    isBurning_k = (snap(k).t <= burn.t_burn + 1e-12);
    if isBurning_k
        Tfill = BC.T_gas; % this should match BC.T_gas if you pass BC in; keeping constant to mimic original
    else
        Tfill = BC.T_amb;  % match BC.T_amb similarly
    end

    burnedMask = r_ref < r_min_k - 1e-12;
    TbarMat(k, burnedMask) = Tfill;
    TmaxMat(k, burnedMask) = Tfill;

    solidMask = ~burnedMask;
    if any(solidMask)
        TbarMat(k, solidMask) = interp1(r_k, Tbar_k, r_ref(solidMask), 'linear', 'extrap');
        TmaxMat(k, solidMask) = interp1(r_k, Tmax_k, r_ref(solidMask), 'linear', 'extrap');
    end
end

% ---- 3D surface plot: Temperature History (radius vs time) ----
[X, Y] = meshgrid(r_rel, t_vec);

figure;
surf(X, Y, TbarMat, 'FaceColor','interp', 'EdgeColor','none');
hold on;

% coarse overlay grid (similar vibe to original)
skipR = max(1, round(Nr_ref/80));
skipT = max(1, round(Nt/60));
mesh(X(1:skipT:end, 1:skipR:end), Y(1:skipT:end, 1:skipR:end), TbarMat(1:skipT:end, 1:skipR:end), ...
    'EdgeColor',[0.6 0.6 0.6], 'EdgeAlpha',0.30, 'FaceColor','none');

% translucent red planes at material interfaces (same as original)
prop_thick  = layers(1).thickness;
liner_thick = layers(2).thickness;
r_if_prop_liner = prop_thick;
r_if_liner_case = prop_thick + liner_thick;
interfaces = [r_if_prop_liner, r_if_liner_case];

Tmin = min(TbarMat(:));
Tmax = max(TbarMat(:));
t_min = min(t_vec);
t_max = max(t_vec);

for kk = 1:numel(interfaces)
    r_int = interfaces(kk);      % already measured from inner surface
    Xp = r_int * ones(2,2);
    Yp = [t_min t_max; t_min t_max];
    Zp = [Tmin Tmin; Tmax Tmax];
    surf(Xp, Yp, Zp, ...
        'FaceColor', [1 0 0], 'FaceAlpha', 0.12, ...
        'EdgeColor', [1 0 0], 'EdgeAlpha', 0.4, 'LineWidth', 1.0);
end

hold off;
colormap(turbo);
colorbar;
grid on; box on;
xlabel('Radius from Original Inner Surface [m]');
ylabel('Time [s]');
zlabel('Temperature [K]');
title('Temperature History (theta-averaged; burned region filled like 1D)');
set(gca, 'XDir','reverse', 'YDir','reverse');

% ---- 2D contour: Temperature Evolution (radius vs time) ----
figure;
contourf(X, Y, TbarMat, 30, 'LineColor','none');
colorbar;
xlabel('Radius from Original Inner Surface [m]');
ylabel('Time [s]');
title('Temperature Evolution (theta-averaged; radius vs time)');
set(gca, 'XDir','reverse', 'YDir','reverse');

% ---- Multi-radii time traces (like original) ----
figure;
hold on;

% Outer surface
plot(t_vec, TbarMat(:, end), 'LineWidth', 2.0, 'DisplayName', 'Outer surface (avg)');

% A few interior radii
r_sample_rel = [0.25, 0.50, 0.75] * r_rel(end);
for kk = 1:numel(r_sample_rel)
    [~, idxr] = min(abs(r_rel - r_sample_rel(kk)));
    plot(t_vec, TbarMat(:, idxr), 'LineWidth', 1.6, ...
        'DisplayName', sprintf('r = %.3f m (avg)', r_rel(idxr)));
end

hold off;
grid on;
xlabel('Time [s]');
ylabel('Temperature [K]');
title('Temperature vs Time at Selected Radii (theta-averaged)');
legend('Location', 'best');

% ---- 2D cross-section plot at final time (your existing 2D view) ----
S = snap(end);
Tlast = S.T;
rlast = S.r_centers(:);
thlast = S.theta(:);

[Rm, THm] = ndgrid(rlast, thlast);
[Xm, Ym] = pol2cart(THm, Rm);

figure;
contourf(Xm, Ym, Tlast, 30, 'LineColor','none');
axis equal tight;
colorbar;
xlabel('x [m]'); ylabel('y [m]');
title(sprintf('Temperature cross-section at t = %.2f s', S.t));

Tin_avg = zeros(Nt,1);
for k = 1:Nt
    Tin = snap(k).T(1,:);
    Tin_avg(k) = mean(Tin);
end
figure;
plot(t_vec, Tin_avg, 'LineWidth', 2);
grid on;
xlabel('Time [s]'); ylabel('Inner solid surface T [K]');
title('Inner solid surface temperature (avg around \theta)');


end

function E = solid_internal_energy(T, r_faces, dtheta, dz, rhoCp_ring)

% Supports T as Nr x Ntheta (2D) or Nr x Ntheta x Nz (3D). dz is axial cell length.
% Total internal energy proxy: sum(rho*cp*T*V) [J], per unit length L
% T is Nr x Ntheta, rhoCp_ring is Nr x 1
Nr = size(T,1);
%Ntheta = size(T,2);

E = 0.0;
for i = 1:Nr
    r_imh = r_faces(i);
    r_iph = r_faces(i+1);
    Vring = 0.5*(r_iph^2 - r_imh^2) * dtheta * dz;   % volume per cell in theta
    % sum over theta for ring i
    E = E + rhoCp_ring(i) * Vring * sum(T(i,:,:),'all');
end
end

% function [dQ_in, dQ_out] = boundary_heat_step(Tprev, Tnew, r_faces, dtheta, dz, BC, T_gas, dt)
% % Compute boundary heat over one timestep using the SAME convection + linearized radiation model
% % as the solver (radiation linearized about Tprev wall temperature).
% %
% % Inputs:
% %   Tprev, Tnew : Nr x Ntheta (or Nr x Ntheta x Nz)
% % Returns:
% %   dQ_in  : positive INTO solid at inner radius
% %   dQ_out : positive LEAVING solid at outer radius
%
% Nr = size(Tprev,1);
% Ntheta = size(Tprev,2);
% Nz = max(1, size(Tprev,3));
%
% dQ_in  = 0.0;
% dQ_out = 0.0;
%
% Ain_cell  = r_faces(1)   * dtheta * dz;
% Aout_cell = r_faces(end) * dtheta * dz;
%
% for kk = 1:Nz
%     for j = 1:Ntheta
%         % inner
%         Tw0 = Tprev(1,j,kk);
%         [h_eff, q_const] = rad_lin_coeff(BC.eps_in, BC.sigma, Tw0, T_gas);
%         h_total = BC.h_in + h_eff;
%         dQ_in = dQ_in + dt * ( h_total*Ain_cell*(T_gas - Tnew(1,j,kk)) + q_const*Ain_cell );
%
%         % outer
%         Tw0 = Tprev(Nr,j,kk);
%         [h_eff, q_const] = rad_lin_coeff(BC.eps_out, BC.sigma, Tw0, BC.T_amb);
%         h_total = BC.h_out + h_eff;
%         dQ_into_solid_outer = dt * ( h_total*Aout_cell*(BC.T_amb - Tnew(Nr,j,kk)) + q_const*Aout_cell );
%         dQ_out = dQ_out - dQ_into_solid_outer;
%     end
% end
% end

% function [dQ_in, dQ_out] = boundary_heat_from_Ab(Tnew, Tprev, r_faces, dtheta, dz, BC, T_gas, dt)
% % Compute boundary heat using the SAME linearized coefficients as assembly,
% % evaluated consistently with the BE solve at Tnew.
% %
% % dQ_in  = positive heat INTO solid at inner boundary
% % dQ_out = positive heat LEAVING solid at outer boundary
%
% Nr = size(Tnew,1);
% Ntheta = size(Tnew,2);
%
% r_in  = r_faces(1);
% r_out = r_faces(end);
%
% Ain  = r_in  * dtheta * dz;
% Aout = r_out * dtheta * dz;
%
% dQ_in  = 0.0;
% dQ_out = 0.0;
%
% % inner boundary (environment = T_gas)
% for j = 1:Ntheta
%     Tw0 = Tprev(1,j);                % linearization point (matches assembly)
%     [h_rad, q_const] = rad_lin_coeff(BC.eps_in, BC.sigma, Tw0, T_gas);
%     h_tot = BC.h_in + h_rad;
%
%     K = h_tot * Ain;                 % diagonal contribution
%     F = h_tot * Ain * T_gas + q_const * Ain;   % RHS contribution
%
%     q_in = (F - K*Tnew(1,j));         % W (since K is W/K, F is W)
%     dQ_in = dQ_in + q_in * dt;        % J
% end
%
% % outer boundary (environment = Tamb)
% for j = 1:Ntheta
%     Tw0 = Tprev(Nr,j);
%     [h_rad, q_const] = rad_lin_coeff(BC.eps_out, BC.sigma, Tw0, BC.T_amb);
%     h_tot = BC.h_out + h_rad;
%
%     K = h_tot * Aout;
%     F = h_tot * Aout * BC.T_amb + q_const * Aout;
%
%     q_into_solid = (F - K*Tnew(Nr,j));    % W into solid (negative when cooling)
%     q_out = max(0, -q_into_solid);        % W leaving solid
%     dQ_out = dQ_out + q_out * dt;
% end
% end

function plot_energy_diagnostics(time, closure_hist, relerr_hist, willDel_hist, Qcond12A_hist, cfg)
%PLOT_ENERGY_DIAGNOSTICS Single-figure diagnostics (tiled subplots).
%
% Panels (always):
%   1) Relative closure error
%   2) Signed residual R
%   3) Absolute residual |R|
%   4) Cumulative residual (drift)
%   5) |R| on log scale
%
% Optional deletion diagnostics (if cfg.doDeletionDiagnostics):
%   6) Deletion-step scatter: R vs Qcond12
%   7) Residual collapse: R and R* = R + alpha*Qcond12
%   8) (optional) Time trace of Qcond12 (commented by default)

dt = time.dt;

% Robust time vector matching residual length
n = numel(closure_hist);
t_step = (0:n-1)' * dt;

R = closure_hist(:);
relerr = relerr_hist(:);

% Finite versions for stats
R_finite = R(isfinite(R));
rel_finite = relerr(isfinite(relerr));
if isempty(R_finite); R_finite = 0; end
if isempty(rel_finite); rel_finite = 0; end

fprintf('Energy closure: RMS(|dE-dQ|)=%.3g J, max relerr=%.3g\n', ...
    rms(R_finite), max(abs(rel_finite)));

% Prep residual with NaNs removed for operations where needed
R0 = R;
R0(~isfinite(R0)) = 0;
cumR = cumsum(R0);

% Determine if deletion diagnostics are available and enabled
doDel = isfield(cfg,'doDeletionDiagnostics') && cfg.doDeletionDiagnostics;

haveQ = ~isempty(Qcond12A_hist) && numel(Qcond12A_hist) == n;
haveW = ~isempty(willDel_hist)  && numel(willDel_hist)  == n;

doDel = doDel && haveQ && haveW;

if isfield(cfg,'doDeletionDiagnostics') && cfg.doDeletionDiagnostics && ~doDel
    fprintf('Deletion diagnostics requested but skipped (missing/mismatched Qcond12A_hist or willDel_hist).\n');
end

% Choose layout
if doDel
    baseTitle = 'Energy diagnostics (with deletion diagnostics)';
else
    baseTitle = 'Energy diagnostics';
end

if isfield(cfg,'diagLabel') && ~isempty(cfg.diagLabel)
    figTitle = sprintf('%s | %s', baseTitle, cfg.diagLabel);
else
    figTitle = baseTitle;
end


fig = figure('Name', figTitle);
clf(fig);

if doDel
    tl = tiledlayout(fig, 3, 3, 'TileSpacing','compact', 'Padding','compact');
else
    tl = tiledlayout(fig, 2, 3, 'TileSpacing','compact', 'Padding','compact');
end

% --- 1) Relative closure error ---
nexttile(tl);
plot(t_step, relerr, 'LineWidth', 1.4);
grid on;
xlabel('Time [s]');
ylabel('Rel err');
title('Relative closure error');

% --- 2) Signed residual ---
nexttile(tl);
plot(t_step, R, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('R [J/step]');
title('Signed residual');

% --- 3) Absolute residual ---
nexttile(tl);
plot(t_step, abs(R), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('|R| [J/step]');
title('Absolute residual');

% --- 4) Cumulative residual ---
nexttile(tl);
plot(t_step, cumR, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('\Sigma R [J]');
title('Cumulative drift');

% --- 5) Absolute residual (log) ---
nexttile(tl);
semilogy(t_step, max(abs(R0), 1e-18), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('|R| [J/step]');
title('|R| (log)');

% If no deletion diagnostics, add a compact text summary panel
if ~doDel
    nexttile(tl);
    axis off;
    txt = sprintf(['RMS(|dE-dQ|) = %.3g J\n' ...
        'max |relerr| = %.3g\n' ...
        'Nsteps = %d\n' ...
        'dt = %.4g s'], ...
        rms(R_finite), max(abs(rel_finite)), n, dt);
    text(0, 0.9, txt, 'FontName','Consolas', 'FontSize', 10, 'VerticalAlignment','top');
end

% --- Deletion diagnostics ---
if doDel
    Q = Qcond12A_hist(:);
    willDel = willDel_hist(:);

    idx = willDel & isfinite(R) & isfinite(Q) & abs(Q) > 1e-12;

    if any(idx)
        ratio = R(idx) ./ (Q(idx) + 1e-12);
        C = corrcoef(R(idx), Q(idx));
        alpha = -mean(R(idx) ./ Q(idx));

        fprintf('Deletion-step match: median R/Qcond12 = %.3g, mean = %.3g\n', ...
            median(ratio), mean(ratio));
        fprintf('corrcoef(R,Qcond12) = %.3f\n', C(1,2));
    else
        alpha = 0.0;
        fprintf('Deletion-step diagnostics: no valid deletion steps found for stats.\n');
    end

    Rstar = R + alpha * Q;

    % --- 6) Scatter: R vs Qcond12 on deletion steps ---
    nexttile(tl);
    if any(idx)
        plot(Q(idx), R(idx), '.', 'MarkerSize', 8);
        grid on;
        xlabel('Qcond12 [J/step]');
        ylabel('R [J/step]');
        title('Deletion steps: R vs Qcond12');
    else
        axis off;
        text(0, 0.6, 'No valid deletion-step points for scatter.', ...
            'FontName','Consolas', 'FontSize', 10);
    end

    % --- 7) Residual collapse test ---
    nexttile(tl);
    plot(t_step, R, 'LineWidth', 1.2); hold on;
    plot(t_step, Rstar, 'LineWidth', 1.2);
    grid on;
    xlabel('Time [s]');
    ylabel('Residual [J/step]');
    title(sprintf('Residual collapse (\\alpha=%.3g)', alpha));
    legend('R','R*', 'Location','best');
    hold off;

    % --- 8) Optional: time trace of Qcond12 ---
    nexttile(tl);
    plot(t_step, Q, 'LineWidth', 1.2);
    grid on;
    xlabel('Time [s]');
    ylabel('Qcond12 [J/step]');
    title('Qcond12 time trace');

    % --- 9) Text summary panel ---
    nexttile(tl);
    axis off;
    if any(idx)
        txt = sprintf(['RMS(|dE-dQ|) = %.3g J\n' ...
            'max |relerr| = %.3g\n' ...
            'corr(R,Qcond12) = %.3f\n' ...
            'median R/Q = %.3g\n' ...
            'alpha = %.3g\n' ...
            'Nsteps = %d, dt=%.4g s'], ...
            rms(R_finite), max(abs(rel_finite)), C(1,2), median(ratio), alpha, n, dt);
    else
        txt = sprintf(['RMS(|dE-dQ|) = %.3g J\n' ...
            'max |relerr| = %.3g\n' ...
            'alpha = %.3g\n' ...
            'Nsteps = %d, dt=%.4g s'], ...
            rms(R_finite), max(abs(rel_finite)), alpha, n, dt);
    end
    text(0, 0.95, txt, 'FontName','Consolas', 'FontSize', 10, 'VerticalAlignment','top');
end

% Global title
title(tl, figTitle);

end

function [h_eff, q_const] = rad_lin_coeff(eps, sigma, Tw0, Tenv)
%RAD_LIN_COEFF Linearized radiation: q_rad = eps*sigma*(Tenv^4 - Tw^4)
% Linearize about Tw0:
%   q_rad ≈ h_rad*(Tenv - Tw) + q_const
% where
%   h_rad   = 4*eps*sigma*Tw0^3
%   q_const = eps*sigma*(Tenv^4 - Tw0^4) - h_rad*(Tenv - Tw0)

    Tw0 = max(Tw0, 1.0);  % avoid 0^3 issues; supports scalar or vector
    Tenv4 = Tenv.^4;      % scalar or vector safe
    Tw04  = Tw0.^4;

    h_rad = 4 * eps * sigma .* (Tw0.^3);

    % constant term to keep linearization exact at Tw=Tw0
    q_const = eps * sigma .* (Tenv4 - Tw04) - h_rad .* (Tenv - Tw0);

    h_eff = h_rad;
end

function hm = harmonic_mean(a, b)
%HARMONIC_MEAN Safe harmonic mean (handles zeros)
den = a + b;
if den <= 0
    hm = 0;
else
    hm = 2 * a * b / den;
end
end

function Ecap = endcap_internal_energy(Tcap0, TcapL, r_faces, cfg, BC)
Ecap = 0.0;

if ~(isfield(cfg,'endcap') && isfield(cfg.endcap,'enable') && cfg.endcap.enable)
    return;
end

% Match assembler geometry exactly
r_inner = r_faces(1);
r_outer = r_faces(end);
A_annulus = pi * (r_outer^2 - r_inner^2);

rhoCp_cap = cfg.endcap.rho * cfg.endcap.cp;
Vcap = A_annulus * cfg.endcap.thickness;   % EXACTLY as in assembly

if ~isfield(cfg.endcap,'ends'); cfg.endcap.ends = 'both'; end
if strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'z0')
    Ecap = Ecap + rhoCp_cap * Vcap * Tcap0;
end
if strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'zL')
    Ecap = Ecap + rhoCp_cap * Vcap * TcapL;
end
end

function Mdiag = build_Mdiag_wall_endcap(r_faces, r_centers, dtheta, dz, rhoCp_ring, dt, cfg)
% Returns Mdiag (J/K) for all unknowns in the SAME ordering as Tvec:
% first Nwall entries are wall cells (sub2ind ordering), then optional endcap nodes.

Nr = numel(r_centers);
Ntheta = cfg.Ntheta;
Nz = cfg.Nz;

Nwall = Nr*Ntheta*Nz;

Mdiag = zeros(Nwall,1);

% Wall cells: M = rhoCp * Vcv (J/K)
idxp = 0;
for kk = 1:Nz
    for j = 1:Ntheta
        for i = 1:Nr
            idxp = idxp + 1;
            r_imh = r_faces(i);
            r_iph = r_faces(i+1);
            Vcv = 0.5*(r_iph^2 - r_imh^2) * dtheta * dz;
            Mdiag(idxp) = rhoCp_ring(i) * Vcv;
        end
    end
end

% Endcaps: lumped M = rhoCp * Vcap (J/K), appended after wall unknowns
useEndcap = isfield(cfg,'endcap') && isfield(cfg.endcap,'enable') && cfg.endcap.enable;
if useEndcap
    if ~isfield(cfg.endcap,'ends'); cfg.endcap.ends = 'both'; end

    r_inner = r_faces(1);
    r_outer = r_faces(end);
    A_annulus = pi*(r_outer^2 - r_inner^2);

    rhoCp_cap = cfg.endcap.rho * cfg.endcap.cp;
    Vcap = A_annulus * cfg.endcap.thickness;
    Mcap = rhoCp_cap * Vcap;

    cap_z0 = strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'z0');
    cap_zL = strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'zL');

    if cap_z0 && cap_zL
        Mdiag = [Mdiag; Mcap; Mcap];
    elseif cap_z0 || cap_zL
        Mdiag = [Mdiag; Mcap];
    end
end

end

function diag = energy_diag_from_Ab(A, b, Tnew_vec, Tprev_vec, Mdiag, dt, Nwall)
% Computes the net heat implied by the assembled linear system:
%   dQ_implied = dt * sum( beff - Aeff*Tnew )
% where Aeff = A - diag(M/dt), beff = b - (M/dt)*Tprev.
%
% Also splits wall vs endcap contributions by row-summing over those equations.

N = numel(Tnew_vec);
assert(numel(Mdiag) == N, 'Mdiag length must match Tvec length.');

Aeff = A - spdiags(Mdiag/dt, 0, N, N);
beff = b - (Mdiag/dt).*Tprev_vec;

rhs_minus_lhs = beff - (Aeff*Tnew_vec);   % W (net into each equation node)
dQ_implied_total = dt * sum(rhs_minus_lhs);  % J

% Energy change from masses (should match dQ_implied_total if BE assembly is self-consistent)
dE_mass_total = sum(Mdiag .* (Tnew_vec - Tprev_vec)); % J

% Split by equation groups (row sets)
rows_wall = 1:Nwall;
dQ_implied_wall = dt * sum(rhs_minus_lhs(rows_wall));
dE_mass_wall    = sum(Mdiag(rows_wall) .* (Tnew_vec(rows_wall) - Tprev_vec(rows_wall)));

rows_cap = (Nwall+1):N;
if isempty(rows_cap)
    dQ_implied_cap = 0;
    dE_mass_cap    = 0;
else
    dQ_implied_cap = dt * sum(rhs_minus_lhs(rows_cap));
    dE_mass_cap    = sum(Mdiag(rows_cap) .* (Tnew_vec(rows_cap) - Tprev_vec(rows_cap)));
end

diag.dQ_implied_total = dQ_implied_total;
diag.dE_mass_total    = dE_mass_total;
diag.resid_total      = dE_mass_total - dQ_implied_total;

diag.dQ_implied_wall  = dQ_implied_wall;
diag.dE_mass_wall     = dE_mass_wall;
diag.resid_wall       = dE_mass_wall - dQ_implied_wall;

diag.dQ_implied_cap   = dQ_implied_cap;
diag.dE_mass_cap      = dE_mass_cap;
diag.resid_cap        = dE_mass_cap - dQ_implied_cap;

end

function cache = build_BE_polar3D_cache(r_faces, r_centers, dtheta, dz, k_ring, rhoCp_ring, dt, cfg, Lg, Ntheta, Nz)
Nr = numel(r_centers);
Nwall = Nr*Ntheta*Nz;
p_all = reshape(1:Nwall, [Nr, Ntheta, Nz]);  % linear indices

useEndcap = isfield(cfg,'endcap') && isfield(cfg.endcap,'enable') && cfg.endcap.enable;
cap_z0 = false; cap_zL = false;
if useEndcap
    if ~isfield(cfg.endcap,'ends'); cfg.endcap.ends = 'both'; end
    cap_z0 = strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'z0');
    cap_zL = strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'zL');
end
Ncap = double(cap_z0) + double(cap_zL);
N = Nwall + Ncap;

% Endcap node indices
p_cap0 = []; p_capL = [];
if cap_z0, p_cap0 = Nwall + 1; end
if cap_zL, p_capL = Nwall + 1 + double(cap_z0); end

% Precompute z centers (used for insul masks and grain-present mask later)
z_centers = ((1:Nz) - 0.5) * dz;

% Surface indexing helpers (wall nodes only)
% All (j,kk) for i=1 and i=Nr
[jGrid,kGrid] = ndgrid(1:Ntheta, 1:Nz);
jList = jGrid(:);
kList = kGrid(:);

idx_i1  = p_all(1,:,:);  idx_i1  = idx_i1(:);
idx_iNr = p_all(Nr,:,:); idx_iNr = idx_iNr(:);
cache.kList_i1 = kList;

% ---- aP_wall: transient diagonal per wall node (vectorized, correct shape) ----
r_imh = r_faces(1:Nr);
r_iph = r_faces(2:Nr+1);

r_imh = r_imh(:);     % force [Nr x 1]
r_iph = r_iph(:);     % force [Nr x 1]

Vcv_i = 0.5*(r_iph.^2 - r_imh.^2) * dtheta * dz;   % [Nr x 1]

rhoCp_ring = rhoCp_ring(:);
k_ring     = k_ring(:);
assert(numel(rhoCp_ring) == Nr, 'rhoCp_ring must be length Nr. Got %d, Nr=%d', numel(rhoCp_ring), Nr);
assert(numel(k_ring)     == Nr, 'k_ring must be length Nr. Got %d, Nr=%d', numel(k_ring), Nr);
aP_i  = (rhoCp_ring(:) .* Vcv_i) / dt;                  % [Nr x 1]
aP_wall_3D = repmat(aP_i(:), 1, Ntheta, Nz);   % [Nr x Ntheta x Nz]
aP_wall    = aP_wall_3D(:);                   % [Nwall x 1]
assert(numel(aP_wall) == Nwall, 'aP_wall numel mismatch: %d vs Nwall %d', numel(aP_wall), Nwall);


% Build A_base (conduction + transient + endcap coupling only)
% Estimate nonzeros per row ~ 10..14
nz_est = N * 14;
I = zeros(nz_est,1);
J = zeros(nz_est,1);
V = zeros(nz_est,1);
idx = 0;

% theta periodic
jp = [2:Ntheta 1];
jm = [Ntheta 1:Ntheta-1];

for kk=1:Nz
    for j=1:Ntheta
        for i=1:Nr
            p = p_all(i,j,kk);

            r_imh = r_faces(i);
            r_iph = r_faces(i+1);
            ri    = r_centers(i);

            % -------- Radial conduction interior only (no BC terms) --------
            if i > 1
                r_face = r_imh;
                Aface  = r_face * dtheta * dz;
                dW     = ri - r_centers(i-1);
                k_face = harmonic_mean(k_ring(i-1), k_ring(i));
                GW     = k_face * Aface / dW;

                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx)+GW;
                q = p_all(i-1,j,kk);
                idx=idx+1; I(idx)=p; J(idx)=q; V(idx)=V(idx)-GW;
            end
            if i < Nr
                r_face = r_iph;
                Aface  = r_face * dtheta * dz;
                dE     = r_centers(i+1) - ri;
                k_face = harmonic_mean(k_ring(i), k_ring(i+1));
                GE     = k_face * Aface / dE;

                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx)+GE;
                q = p_all(i+1, j, kk);
                idx=idx+1; I(idx)=p; J(idx)=q; V(idx)=V(idx)-GE;
            end

            % -------- Theta conduction (periodic) --------
            dr_i   = (r_iph - r_imh);
            Atheta = dr_i * dz;
            ds     = ri * dtheta;
            Gth    = k_ring(i) * Atheta / ds;

            jE = jp(j); jW = jm(j);
            qE = p_all(i, jE, kk);
            qW = p_all(i, jW, kk);

            idx=idx+1; I(idx)=p; J(idx)=p;  V(idx)=V(idx) + 2*Gth;
            idx=idx+1; I(idx)=p; J(idx)=qE; V(idx)=V(idx) - Gth;
            idx=idx+1; I(idx)=p; J(idx)=qW; V(idx)=V(idx) - Gth;

            % -------- Axial conduction interior; endcap coupling if enabled --------
            Az = 0.5*(r_iph^2 - r_imh^2) * dtheta;
            Gz = k_ring(i) * Az / dz;

            if kk > 1
                q = p_all(i,j,kk-1);
                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx)+Gz;
                idx=idx+1; I(idx)=p; J(idx)=q; V(idx)=V(idx)-Gz;
            else
                if cap_z0
                    if ~isfield(cfg.endcap,'thickness'); cfg.endcap.thickness = dz; end
                    if ~isfield(cfg.endcap,'k'); cfg.endcap.k = 1.0; end
                    Gcap = cfg.endcap.k * Az / cfg.endcap.thickness;

                    idx=idx+1; I(idx)=p;      J(idx)=p;      V(idx)=V(idx) + Gcap;
                    idx=idx+1; I(idx)=p;      J(idx)=p_cap0; V(idx)=V(idx) - Gcap;

                    idx=idx+1; I(idx)=p_cap0; J(idx)=p_cap0; V(idx)=V(idx) + Gcap;
                    idx=idx+1; I(idx)=p_cap0; J(idx)=p;      V(idx)=V(idx) - Gcap;
                end
            end

            if kk < Nz
                q = p_all(i, j, kk+1);
                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx)+Gz;
                idx=idx+1; I(idx)=p; J(idx)=q; V(idx)=V(idx)-Gz;
            else
                if cap_zL
                    if ~isfield(cfg.endcap,'thickness'); cfg.endcap.thickness = dz; end
                    if ~isfield(cfg.endcap,'k'); cfg.endcap.k = 1.0; end
                    Gcap = cfg.endcap.k * Az / cfg.endcap.thickness;

                    idx=idx+1; I(idx)=p;      J(idx)=p;      V(idx)=V(idx) + Gcap;
                    idx=idx+1; I(idx)=p;      J(idx)=p_capL; V(idx)=V(idx) - Gcap;

                    idx=idx+1; I(idx)=p_capL; J(idx)=p_capL; V(idx)=V(idx) + Gcap;
                    idx=idx+1; I(idx)=p_capL; J(idx)=p;      V(idx)=V(idx) - Gcap;
                end
            end

            % -------- Transient diagonal --------
            idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + aP_wall(p);
        end
    end
end

% Endcap transient diagonal gets included in A_base too
aPcap = 0;
A_annulus = [];
if useEndcap
    if ~isfield(cfg.endcap,'thickness'); cfg.endcap.thickness = dz; end
    if ~isfield(cfg.endcap,'rho'); cfg.endcap.rho = 1000; end
    if ~isfield(cfg.endcap,'cp');  cfg.endcap.cp  = 1000; end

    r_inner = r_faces(1);
    r_outer = r_faces(end);
    A_annulus = pi * (r_outer^2 - r_inner^2);
    Vcap = A_annulus * cfg.endcap.thickness;
    aPcap = (cfg.endcap.rho * cfg.endcap.cp) * Vcap / dt;

    if cap_z0
        idx=idx+1; I(idx)=p_cap0; J(idx)=p_cap0; V(idx)=V(idx) + aPcap;
    end
    if cap_zL
        idx=idx+1; I(idx)=p_capL; J(idx)=p_capL; V(idx)=V(idx) + aPcap;
    end
end

I = I(1:idx); J = J(1:idx); V = V(1:idx);
A_base = sparse(I,J,V,N,N);

% Pack cache
cache = struct();
cache.Nr = Nr; cache.Ntheta = Ntheta; cache.Nz = Nz;
cache.Nwall = Nwall; cache.Ncap = Ncap; cache.N = N;
cache.cap_z0 = cap_z0; cache.cap_zL = cap_zL;
cache.p_cap0 = p_cap0; cache.p_capL = p_capL;

cache.r_faces = r_faces;
cache.r_centers = r_centers;
cache.dtheta = dtheta; cache.dz = dz;
cache.z_centers = z_centers;

cache.A_base = A_base;
assert(isvector(aP_wall), 'aP_wall is not a vector, size=%s', mat2str(size(aP_wall)));
assert(numel(aP_wall) == Nwall, 'aP_wall numel mismatch: %d vs Nwall %d', numel(aP_wall), Nwall);
cache.aP_wall = full(aP_wall(:));

cache.idx_i1 = idx_i1;
cache.idx_iNr = idx_iNr;

cache.Ain_face = r_faces(1) * dtheta * dz;
cache.Aout_face = r_faces(end) * dtheta * dz;

cache.A_annulus = A_annulus;
cache.aPcap = aPcap;
end

function [A, b] = assemble_from_cache_BE_polar3D(cache, Tprev, Tcap0_prev, TcapL_prev, BC, T_gas, cfg, Lg)
% --- Start from cached base ---
A = cache.A_base;

% --- Accumulate ONLY diagonal BC contributions here (size N x 1) ---
diag_add = zeros(cache.N,1);
b_add = zeros(cache.N,1);

N = cache.N;
Nwall = cache.Nwall;

% Start with transient RHS (vectorized)
b = zeros(N,1);
ap = cache.aP_wall(:);
Tp = Tprev(:);
assert(numel(ap) == numel(Tp), 'aP_wall (%d) and Tprev (%d) size mismatch', numel(ap), numel(Tp));
b(1:Nwall) = ap .* Tp;

% Endcap transient RHS
useEndcap = cache.Ncap > 0;
if useEndcap
    if cache.cap_z0
        b(cache.p_cap0) = b(cache.p_cap0) + cache.aPcap * Tcap0_prev;
    end
    if cache.cap_zL
        b(cache.p_capL) = b(cache.p_capL) + cache.aPcap * TcapL_prev;
    end
end

Nr = cache.Nr; Ntheta = cache.Ntheta; Nz = cache.Nz;
z_centers = cache.z_centers;

% Grain present mask
grain_present = (z_centers <= Lg + 1e-12);

% Step 2 insulation multiplier along z
m_in = ones(Nz,1);
useWallInsul = isfield(cfg,'wall_insul') && isfield(cfg.wall_insul,'enable') && cfg.wall_insul.enable;
if useWallInsul
    if ~isfield(cfg.wall_insul,'end');    cfg.wall_insul.end = 'z0'; end
    if ~isfield(cfg.wall_insul,'length'); cfg.wall_insul.length = 0.0; end
    if ~isfield(cfg.wall_insul,'factor'); cfg.wall_insul.factor = 0.0; end

    if strcmpi(cfg.wall_insul.end,'z0')
        mask = (z_centers <= cfg.wall_insul.length);
    elseif strcmpi(cfg.wall_insul.end,'zL')
        mask = (z_centers >= (max(z_centers) - cfg.wall_insul.length));
    else
        error('cfg.wall_insul.end must be ''z0'' or ''zL''.');
    end
    m_in(mask) = cfg.wall_insul.factor;
end

% Nozzle proxy scaling (Step 4)
useNozProxy = isfield(cfg,'noz') && isfield(cfg.noz,'enableProxy') && cfg.noz.enableProxy;
r_inner = cache.r_faces(1);
Aport = pi * r_inner^2;
At = [];
if useNozProxy
    r_t = cfg.noz.r_t;
    At = pi * r_t^2;
end

% ---------- Inner wall BC (i=1) ----------
idx_i1 = cache.idx_i1;                    % length Ntheta*Nz
Tw_i1  = Tprev(idx_i1);
Aface  = cache.Ain_face;

% --- k-index list for i==1 nodes (must match idx_i1 ordering) ---
if isfield(cache,'kList_i1') && ~isempty(cache.kList_i1)
    kk_list = cache.kList_i1(:);
else
    [~, kGrid] = ndgrid(1:Ntheta, 1:Nz);
    kk_list = kGrid(:);
end

inGrain = grain_present(kk_list);     % logical, length(idx_i1)

% =====================
% Grain-present region
% =====================
if any(inGrain)
    pG  = idx_i1(inGrain);
    TwG = Tw_i1(inGrain);
    mG  = m_in(kk_list(inGrain));

    [h_effG, q_constG] = rad_lin_coeff(BC.eps_in, BC.sigma, TwG, T_gas);
    htotG = mG .* (BC.h_in + h_effG);

    active = (htotG > 0);
    if any(active)
        pG = pG(active);
        htotG = htotG(active);
        mG = mG(active);
        q_constG = q_constG(active);

        diag_add(pG) = diag_add(pG) + htotG * Aface;
        b(pG)        = b(pG)        + (htotG * Aface) .* T_gas + (mG .* q_constG) * Aface;
    end
end

% =====================
% Aft cavity region
% =====================
if any(~inGrain)
    pC  = idx_i1(~inGrain);
    TwC = Tw_i1(~inGrain);

    if isfield(BC,'T_gas_cav_mode') && strcmpi(BC.T_gas_cav_mode,'ambient')
        Tenv = BC.T_amb;
    else
        Tenv = T_gas;
    end

    eps_cav = BC.eps_in_cav;

    % h_cav may be constant or nozzle-scaled, but it is the same for all cavity faces
    h_cav = cfg.noz.h_cav_base;
    if useNozProxy
        ratio = cfg.noz.CdA * (At / max(Aport, 1e-12));
        h_cav = cfg.noz.h_cav_base * (ratio ^ cfg.noz.alpha);
        h_cav = min(max(h_cav, cfg.noz.h_cav_min), cfg.noz.h_cav_max);
    end

    [h_effC, q_constC] = rad_lin_coeff(eps_cav, BC.sigma, TwC, Tenv);
    htotC = h_cav + h_effC;

    active = (htotC > 0);
    if any(active)
        pC = pC(active);
        htotC = htotC(active);
        q_constC = q_constC(active);

        diag_add(pC) = diag_add(pC) + htotC * Aface;
        b(pC)        = b(pC)        + (htotC * Aface) .* Tenv + (q_constC * Aface);
    end
end


% ---------- Outer wall BC (i=Nr) ----------
idx_iNr    = cache.idx_iNr;
Tw_iNr     = Tprev(idx_iNr);
Aface_out  = cache.Aout_face;

[h_effO, q_constO] = rad_lin_coeff(BC.eps_out, BC.sigma, Tw_iNr, BC.T_amb);
htotO = BC.h_out + h_effO;

active = (htotO > 0);
if any(active)
    pO = idx_iNr(active);
    htotO = htotO(active);
    q_constO = q_constO(active);

    diag_add(pO) = diag_add(pO) + htotO * Aface_out;
    b(pO)        = b(pO)        + (htotO * Aface_out) .* BC.T_amb + (q_constO * Aface_out);
end

% ---------- Endcap lumped BCs ----------
if useEndcap
    if ~isfield(cfg.endcap,'h_in');    cfg.endcap.h_in = 0; end
    if ~isfield(cfg.endcap,'eps_in');  cfg.endcap.eps_in = 0; end
    if ~isfield(cfg.endcap,'h_out');   cfg.endcap.h_out = 0; end
    if ~isfield(cfg.endcap,'eps_out'); cfg.endcap.eps_out = 0; end

    % Defaults for gas exposure flags
    gas_exposed_z0 = true;
    gas_exposed_zL = true;
    if isfield(cfg.endcap,'gas_exposed_z0'); gas_exposed_z0 = cfg.endcap.gas_exposed_z0; end
    if isfield(cfg.endcap,'gas_exposed_zL'); gas_exposed_zL = cfg.endcap.gas_exposed_zL; end

    A_ann = cache.A_annulus;

    if cache.cap_z0
        pcap = cache.p_cap0;
        Tw_cap = Tcap0_prev;

        % inside to gas
        if gas_exposed_z0
            [h_rad_in, q_const_in] = rad_lin_coeff(cfg.endcap.eps_in, BC.sigma, Tw_cap, T_gas);
            htot_in = cfg.endcap.h_in + h_rad_in;
            if htot_in > 0
                diag_add(pcap) = diag_add(pcap) + htot_in*A_ann;
                b(pcap) = b(pcap) + htot_in*A_ann*T_gas + q_const_in*A_ann;
            end
        end

        % outside to ambient
        [h_rad_out, q_const_out] = rad_lin_coeff(cfg.endcap.eps_out, BC.sigma, Tw_cap, BC.T_amb);
        htot_out = cfg.endcap.h_out + h_rad_out;
        if htot_out > 0
            diag_add(pcap) = diag_add(pcap) + htot_out*A_ann;
            b(pcap) = b(pcap) + htot_out*A_ann*BC.T_amb + q_const_out*A_ann;
        end
    end

    if cache.cap_zL
        pcap = cache.p_capL;
        Tw_cap = TcapL_prev;

        if gas_exposed_zL
            [h_rad_in, q_const_in] = rad_lin_coeff(cfg.endcap.eps_in, BC.sigma, Tw_cap, T_gas);
            htot_in = cfg.endcap.h_in + h_rad_in;
            if htot_in > 0
                diag_add(pcap) = diag_add(pcap) + htot_in*A_ann;
                b(pcap) = b(pcap) + htot_in*A_ann*T_gas + q_const_in*A_ann;
            end
        end

        [h_rad_out, q_const_out] = rad_lin_coeff(cfg.endcap.eps_out, BC.sigma, Tw_cap, BC.T_amb);
        htot_out = cfg.endcap.h_out + h_rad_out;
        if htot_out > 0
            diag_add(pcap) = diag_add(pcap) + htot_out*A_ann;
            b(pcap) = b(pcap) + htot_out*A_ann*BC.T_amb + q_const_out*A_ann;
        end
    end
end

% Final system
A = cache.A_base + spdiags(diag_add, 0, N, N);
end

function x = pack_state_vec(Twall, cap_z0, cap_zL, Tcap0, TcapL)
% Pack [Twall(:); optional caps] in the same ordering as assembler expects.
x = Twall(:);
if cap_z0, x = [x; Tcap0]; end
if cap_zL, x = [x; TcapL]; end
end

function [Twall, Tcap0, TcapL] = unpack_state_vec(x, Nr, Ntheta, Nz, cap_z0, cap_zL)
% Unpack x into wall and optional caps.
Nwall = Nr*Ntheta*Nz;
Twall = reshape(x(1:Nwall), [Nr, Ntheta, Nz]);
p = Nwall;

Tcap0 = [];
TcapL = [];
if cap_z0
    Tcap0 = x(p+1); p = p+1;
end
if cap_zL
    TcapL = x(p+1);
end
end

function [isBurning,T_gas,ndel,willDelete,willDel_hist,Lg,grain_present,n_grain] = compute_step_state(step, t, T, r_centers, dz, burn, geom, BC, cfg)

% ---- burn state + gas temp ----
isBurning = (t < burn.t_burn + 1e-12);
if isBurning
    T_gas = BC.T_gas;
else
    tau_gas = 20;  % [s] tune to test data
    T_gas = BC.T_amb + (BC.T_gas - BC.T_amb) * exp(-(t - burn.t_burn)/tau_gas);
end

% ---- Decide deletion NOW (based on current t) ----
if isBurning
    r_front = geom.r_inner0 + burn.rburn * t;
    ndel = sum(r_centers < r_front);
else
    ndel = 0;
end
ndel = min(ndel, numel(r_centers)-2);

willDelete = (ndel > 0);
willDel_hist(step+1) = willDelete;

if mod(step,cfg.status_every)==0 && willDelete
    fprintf('   deletion planned: ndel = %d (Nr before = %d)\n', ndel, numel(r_centers));
end

% ---- current grain length Lg(t) ----
if isfield(cfg,'grain') && isfield(cfg.grain,'enable_axial_regression') && cfg.grain.enable_axial_regression
    Lg = max(0.0, cfg.grain.Lg0 - cfg.grain.rb_end * t);
else
    Lg = cfg.L;
end

Nz = size(T,3);
z_centers = ((1:Nz) - 0.5) * dz;
grain_present = (z_centers <= Lg + 1e-12);
n_grain = sum(grain_present);

end