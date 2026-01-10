function SRM_Heat_2D_Polar()
clc; clear;

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

liner.name = 'EPDM';
liner.k = 0.20;
liner.rho = 1100;
liner.cp = 1500;

caseCF.name = 'CF_composite';
caseCF.k = 1.0;
caseCF.rho = 1600;
caseCF.cp = 900;

% layers
layers = struct([]);
layers(1).material = prop;
layers(1).thickness = 0.0143;
layers(2).material = liner;
layers(2).thickness = 0.002;
layers(3).material = caseCF;
layers(3).thickness = 0.00338;

% geometry / burn
geom.r_inner0 = 0.050;
burn.t_burn = 6;
burn.rburn  = layers(1).thickness / burn.t_burn;  % aligns with your intent

% boundary conditions
BC.T_gas = 2300;          % "inner" gas temperature
BC.T_amb = 300;           % ambient
BC.h_out = 30;            % outer convection coefficient
BC.h_in  = 2000;          % inner convection coefficient (set realistically later)
BC.eps_in  = 0.8;
BC.eps_out = 0.8;
BC.sigma = 5.670374419e-8;

time.dt = 0.010;
time.t_final = 100.0;

cfg.store_every = 10;     % store every N steps
cfg.max_cells = 2e6;      % safety
cfg.Ntheta = 64;          % circumferential resolution (>= 8)
cfg.L = 1.0;              % unit length (m)

% ---- runtime status monitoring ----
cfg.status_every = 20;     % print every N time steps
cfg.time_warn    = 2.0;    % warn if a step takes > this many seconds
cfg.doEnergyDiagnostics = false;     % master switch for all plots/prints below
cfg.doDeletionDiagnostics = false;   % requires Qcond12A_hist + willDel_hist; set false to skip those plots

%% -------------------- Build polar FV grid --------------------
dr = burn.rburn * time.dt;         % radial step aligned to burn per step
[dtheta, theta_centers] = uniform_theta_grid(cfg.Ntheta);

[r_faces, r_centers, layer_idx_of_ring] = build_radial_faces_with_interfaces(geom, layers, dr);
Nr = numel(r_centers);
Ntheta = cfg.Ntheta;

nCells = Nr * Ntheta;
if nCells > cfg.max_cells
    error('Too many cells (%d). Reduce Ntheta or increase dt.', nCells);
end

% material props per radial ring (piecewise constant in theta by default)
[k_ring, rhoCp_ring] = ring_props(layers, layer_idx_of_ring);

% initial condition
T = BC.T_amb * ones(Nr, Ntheta);

% diagnostics: rough Fourier number using smallest ring dr
alpha_prop = prop.k / (prop.rho*prop.cp);
Fo = alpha_prop * time.dt / dr^2;
fprintf('Grid: Nr=%d, Ntheta=%d, dr=%.6g m, dt=%.6g s, Fo_prop~%.3g\n', Nr, Ntheta, dr, time.dt, Fo);

%% -------------------- Time integration --------------------
nsteps = round(time.t_final / time.dt);

% ---- Energy balance tracking ----
E_solid_hist   = nan(nsteps+1,1);   % total internal energy [J] (per unit length if L=1)
dE_hist        = nan(nsteps,1);     % stepwise delta E
dQnet_hist     = nan(nsteps,1);     % stepwise net boundary heat into solid
closure_hist   = nan(nsteps,1);     % dE - dQnet
relerr_hist    = nan(nsteps,1);     % (dE - dQnet)/max(|dQnet|,|dE|)
Qcond12A_hist = nan(nsteps,1);
Erem_hist    = nan(nsteps,1);
willDel_hist = false(nsteps,1);
Qcond12A_hist = nan(nsteps,1);

% store initial energy
E_solid_hist(1) = solid_internal_energy(T, r_faces, dtheta, cfg.L, rhoCp_ring);

times = [];
Tout_avg = [];
Tout_max = [];
snap = struct('t',{},'r_centers',{},'theta',{},'T',{});

for step = 0:nsteps
    t = step * time.dt;

    % ---- store ----
    if mod(step, cfg.store_every)==0
        times(end+1,1) = t;
        [Tavg, Tmax] = outer_wall_stats(T);
        Tout_avg(end+1,1) = Tavg;
        Tout_max(end+1,1) = Tmax;

        snap(end+1).t = t;
        snap(end).r_centers = r_centers;
        snap(end).theta = theta_centers;
        snap(end).T = T;

        % 1D equivalents for original-code style plots
        snap(end).Tbar_r = mean(T, 2);
        snap(end).Tmax_r = max(T, [], 2);
    end

    if step == nsteps
        break;
    end

    % ---- burn state + gas temp ----
    isBurning = (t < burn.t_burn + 1e-12);

    if isBurning
        T_gas = BC.T_gas;
    else
        tau_gas = 20;  % [s] tune to test data
        T_gas = BC.T_amb + (BC.T_gas - BC.T_amb) * exp(-(t - burn.t_burn)/tau_gas);
    end

    % Decide deletion NOW (needed for diagnostics and for later deletion)
    willDelete = isBurning && ~isempty(layer_idx_of_ring) && (layer_idx_of_ring(1) == 1);
    willDel_hist(step+1) = willDelete;


    % ---- assemble + solve ----
    Tprev = T;
    t_start_step = tic;

    t_asm = tic;
    [A, b] = assemble_BE_polar_2D( ...
        Tprev, r_faces, r_centers, dtheta, cfg.L, ...
        k_ring, rhoCp_ring, time.dt, ...
        BC, isBurning, T_gas);
    asm_time = toc(t_asm);

    t_solve = tic;
    Tvec = A \ b;
    solve_time = toc(t_solve);

    T = reshape(Tvec, [Nr, Ntheta]);
    step_time = toc(t_start_step);

    % =========================
    % Boundary heat from (A,b) (single source of truth)
    % =========================
    Vring   = 0.5*(r_faces(2:end).^2 - r_faces(1:end-1).^2) * dtheta * cfg.L; % Nr x 1
    aP_ring = rhoCp_ring(:) .* Vring(:) / time.dt;                           % Nr x 1

    dQ_in  = 0.0;
    dQ_out = 0.0;

    for j = 1:Ntheta
        % ---------- INNER boundary row (i=1) ----------
        i = 1;
        p = (j-1)*Nr + i;

        r_imh = r_faces(i);
        r_iph = r_faces(i+1);
        ri    = r_centers(i);

        % East conductance GE
        r_face = r_iph;
        Aface  = r_face * dtheta * cfg.L;
        dE     = r_centers(2) - ri;
        k_face = harmonic_mean(k_ring(1), k_ring(2));
        GE     = k_face * Aface / dE;

        % Theta conductance Gth
        dr_i   = (r_iph - r_imh);
        Atheta = dr_i * cfg.L;
        ds     = ri * dtheta;
        Gth    = k_ring(1) * Atheta / ds;

        App   = full(A(p,p));
        hA_in = App - (aP_ring(1) + GE + 2*Gth);

        b_bc       = b(p) - aP_ring(1) * Tprev(1,j);   % = hA_in*T_gas + qconstA_in
        qconstA_in = b_bc - hA_in * T_gas;

        dQ_in = dQ_in + time.dt * ( hA_in * (T_gas - T(1,j)) + qconstA_in );

        % ---------- OUTER boundary row (i=Nr) ----------
        i = Nr;
        p = (j-1)*Nr + i;

        r_imh = r_faces(i);
        r_iph = r_faces(i+1);
        ri    = r_centers(i);

        % West conductance GW
        r_face = r_imh;
        Aface  = r_face * dtheta * cfg.L;
        dW     = ri - r_centers(i-1);
        k_face = harmonic_mean(k_ring(i-1), k_ring(i));
        GW     = k_face * Aface / dW;

        % Theta conductance Gth
        dr_i   = (r_iph - r_imh);
        Atheta = dr_i * cfg.L;
        ds     = ri * dtheta;
        Gth    = k_ring(i) * Atheta / ds;

        App    = full(A(p,p));
        hA_out = App - (aP_ring(i) + GW + 2*Gth);

        b_bc        = b(p) - aP_ring(i) * Tprev(i,j);  % = hA_out*T_amb + qconstA_out
        qconstA_out = b_bc - hA_out * BC.T_amb;

        dQ_into_solid_outer = time.dt * ( hA_out * (BC.T_amb - T(i,j)) + qconstA_out );

        % define dQ_out as positive LEAVING solid
        dQ_out = dQ_out - dQ_into_solid_outer;
    end

    dQnet = dQ_in - dQ_out;


    % =========================
    % Definitive boundary heat from (A,b)
    % =========================
    % We decompose each boundary row into:
    %   A_pp = aP + G_radial + 2*Gtheta + (hA)_bc
    % and b_p = aP*Tprev + (hA*Tenv + qconstA)
    % so we can recover (hA)_bc and qconstA exactly and compute:
    %   Q_in_cell  = dt * [ hA*(Tenv - Tnew) + qconstA ]
    %
    % This matches the solver exactly (same linearization, same sign).

    % Geometry + transient coeff per ring
    Vring  = 0.5*(r_faces(2:end).^2 - r_faces(1:end-1).^2) * dtheta * cfg.L; % Nr x 1, per theta cell
    aP_ring = rhoCp_ring(:) .* Vring(:) / time.dt;                           % Nr x 1

    dQin_fromAb  = 0.0;
    dQout_fromAb = 0.0;

    for j = 1:Ntheta
        % ---------- INNER boundary row (i=1) ----------
        i = 1;
        p = (j-1)*Nr + i;

        % Reconstruct GE and Gth exactly like assembler
        r_imh = r_faces(i);
        r_iph = r_faces(i+1);
        ri    = r_centers(i);

        % East neighbor conductance (i=1 to i=2)
        r_face = r_iph;
        Aface  = r_face * dtheta * cfg.L;
        dE     = r_centers(2) - ri;
        k_face = harmonic_mean(k_ring(1), k_ring(2));
        GE     = k_face * Aface / dE;

        % Theta conductance at ring i=1
        dr_i   = (r_iph - r_imh);
        Atheta = dr_i * cfg.L;
        ds     = ri * dtheta;
        Gth    = k_ring(1) * Atheta / ds;

        % Recover hA from diagonal
        App = full(A(p,p));
        hA_in = App - (aP_ring(1) + GE + 2*Gth);

        % Recover boundary source term from b
        b_bc = b(p) - aP_ring(1) * Tprev(1,j);         % = hA*Tgas + qconstA
        qconstA_in = b_bc - hA_in * T_gas;

        % Compute heat INTO solid over dt using T^{n+1} (implicit form)
        dQin_fromAb = dQin_fromAb + time.dt * ( hA_in * (T_gas - T(1,j)) + qconstA_in );

        % ---------- OUTER boundary row (i=Nr) ----------
        i = Nr;
        p = (j-1)*Nr + i;

        r_imh = r_faces(i);
        r_iph = r_faces(i+1);
        ri    = r_centers(i);

        % West neighbor conductance (i=Nr to i=Nr-1)
        r_face = r_imh;
        Aface  = r_face * dtheta * cfg.L;
        dW     = ri - r_centers(i-1);
        k_face = harmonic_mean(k_ring(i-1), k_ring(i));
        GW     = k_face * Aface / dW;

        % Theta conductance at ring i=Nr
        dr_i   = (r_iph - r_imh);
        Atheta = dr_i * cfg.L;
        ds     = ri * dtheta;
        Gth    = k_ring(i) * Atheta / ds;

        App = full(A(p,p));
        hA_out = App - (aP_ring(i) + GW + 2*Gth);

        b_bc = b(p) - aP_ring(i) * Tprev(i,j);         % = hA*Tamb + qconstA
        qconstA_out = b_bc - hA_out * BC.T_amb;

        % This quantity is heat INTO solid from ambient (usually negative when wall>amb)
        dQ_into_solid_outer = time.dt * ( hA_out * (BC.T_amb - T(i,j)) + qconstA_out );

        % For bookkeeping, define dQout as positive LEAVING solid
        dQout_fromAb = dQout_fromAb - dQ_into_solid_outer;
    end

    dQnet_fromAb = dQin_fromAb - dQout_fromAb;

    % ---- Energy balance for this step (BEFORE deleting rings) ----
    % --- Energies computed on the same domain as the solve (pre-deletion) ---
    E_prev = E_solid_hist(step+1);  % energy of the carried state at start of this step
    E_new  = solid_internal_energy(T, r_faces, dtheta, cfg.L, rhoCp_ring);

    % Energy that leaves the computational domain if we delete ring 1 AFTER the solve
    E_removed = 0.0;
    if willDelete
        r_imh = r_faces(1);
        r_iph = r_faces(2);
        Vcell = 0.5*(r_iph^2 - r_imh^2) * dtheta * cfg.L;   % per-theta-cell volume
        E_removed = rhoCp_ring(1) * Vcell * sum(T(1,:));    % sum over theta
    end

    % Corrected energy change over the step for a moving/deleting domain
    dE_corr = (E_new - E_prev);

    % --- IMPORTANT: store the energy of the state that will be carried into the NEXT step ---
    if willDelete
        E_solid_hist(step+2) = E_new - E_removed;   % energy AFTER deletion (matches T after you delete)
    else
        E_solid_hist(step+2) = E_new;
    end

    closure = dE_corr - dQnet;

    % relative error with stable scaling
    Qscale = max(abs(dQ_in) + abs(dQ_out), 1.0);
    relerr = closure / Qscale;

    % store energy histories
    dE_hist(step+1)      = dE_corr;
    dQnet_hist(step+1)   = dQnet;
    closure_hist(step+1) = closure;
    relerr_hist(step+1)  = relerr;

    % ---- NEW diagnostic: Qcond12 computed directly from matrix coupling (optional but useful) ----
    % This is not used to "fix" anything yet; it is for pinpointing.
    Qcond12_fromA = NaN;
    if willDelete && Nr >= 2
        % indices for (i=1,j) and (i=2,j) for all theta columns
        p1 = sub2ind([Nr, Ntheta], ones(1,Ntheta), 1:Ntheta);
        p2 = sub2ind([Nr, Ntheta], 2*ones(1,Ntheta), 1:Ntheta);

        % extract coupling coefficients for each theta: A(p1,p2)
        lin = sub2ind(size(A), p1, p2);
        G12_mat = full(-A(lin)).';          % Ntheta x 1

        dT12 = (T(2,:) - T(1,:)).';         % Ntheta x 1

        Qcond12_fromA = time.dt * sum(G12_mat .* dT12);  % J per step (into ring 1 from ring 2)
    end
    Qcond12A_hist(step+1) = Qcond12_fromA;   % create this array before the loop

    % ---- status output ----
    if mod(step, cfg.status_every) == 0 || step == 1
        Tmax_now = max(T(:));
        Tmin_now = min(T(:));
        Tout_max_now = max(T(end,:));

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
            Nr, Ntheta, Nr*Ntheta, ...
            Tmin_now, Tmax_now, Tout_max_now, ...
            asm_time, solve_time, step_time);

        if step_time > cfg.time_warn
            fprintf('   ⚠ step time exceeded %.1f s\n', cfg.time_warn);
        end
    end

    if mod(step, cfg.status_every) == 0
        fprintf('Progress: %5.1f %%\n', 100*step/nsteps);
    end

    % ---- Move burn front by deleting first radial ring while still in propellant ----
    if willDelete
        T(1,:) = [];
        r_centers(1) = [];
        r_faces(1) = [];
        layer_idx_of_ring(1) = [];

        Nr = numel(r_centers);
        [k_ring, rhoCp_ring] = ring_props(layers, layer_idx_of_ring);
    end

end

plot_results(snap, layers, geom, burn, BC);

if isfield(cfg,'doEnergyDiagnostics') && cfg.doEnergyDiagnostics
    plot_energy_diagnostics(time, closure_hist, relerr_hist, willDel_hist, Qcond12A_hist, cfg);
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
r_outer = r0 + total_thick;

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
layer_idx_of_ring = [];

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

function [A, b] = assemble_BE_polar_2D(Tprev, r_faces, r_centers, dtheta, L, k_ring, rhoCp_ring, dt, BC, isBurning, T_gas)
% Assemble sparse Backward Euler FV system for 2D polar conduction (r,theta).
% Unknown ordering p = (j-1)*Nr + i

Nr = numel(r_centers);
Ntheta = size(Tprev,2);
N = Nr * Ntheta;

% Preallocate sparse with ~5 nonzeros per row (r- and theta-neighbors + diag)
nz_est = N * 12;
I = zeros(nz_est,1);
J = zeros(nz_est,1);
V = zeros(nz_est,1);
b = zeros(N,1);

idx = 0;

% Convenience: periodic theta indexing
jp = @(j) (j < Ntheta) * (j+1) + (j==Ntheta) * 1;
jm = @(j) (j > 1)     * (j-1) + (j==1)     * Ntheta;

% Linearize radiation about previous wall temperature for each boundary cell:
% q_rad = eps*sigma*(T_env^4 - T_w^4) ≈ h_rad*(T_env - T_w) + const
% with h_rad = 4*eps*sigma*Tprev_w^3 and const = eps*sigma*(T_env^4 - Tprev_w^4) - h_rad*(T_env - Tprev_w)
% This gives a consistent "implicit" term in T_w and a source term in b.
%
% For convection + radiation combined: q = h_c*(T_env - T_w) + h_rad*(T_env - T_w) + const
% => q = h_eff*(T_env - T_w) + const, where h_eff = h_c + h_rad

for j = 1:Ntheta
    for i = 1:Nr
        p = (j-1)*Nr + i;

        % Control volume geometry
        r_imh = r_faces(i);       % r_{i-1/2}
        r_iph = r_faces(i+1);     % r_{i+1/2}
        ri    = r_centers(i);

        Vcv = 0.5*(r_iph^2 - r_imh^2) * dtheta * L;  % volume per unit length L
        aP = rhoCp_ring(i) * Vcv / dt;
        bP = aP * Tprev(i,j);

        % ---------- Radial conduction neighbors ----------
        % West (inner) radial face between i-1 and i
        if i > 1
            % face at r_imh, conductance G = k_face * A / d
            r_face = r_imh;
            Aface = r_face * dtheta * L;
            % distances between centers
            dW = ri - r_centers(i-1);
            k_face = harmonic_mean(k_ring(i-1), k_ring(i));
            GW = k_face * Aface / dW;

            % add to matrix
            idx = idx + 1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + GW;
            q = (j-1)*Nr + (i-1);
            idx = idx + 1; I(idx)=p; J(idx)=q; V(idx)=V(idx) - GW;
        else
            % Inner boundary at r = r_faces(1)
            % Robin BC to gas, active during and after burn
            Tw0 = Tprev(i,j);
            [h_eff, q_const] = rad_lin_coeff(BC.eps_in, BC.sigma, Tw0, T_gas);
            h_total = BC.h_in + h_eff;

            r_face = r_imh;
            Aface  = r_face * dtheta * L;

            idx = idx + 1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + h_total*Aface;
            bP = bP + h_total*Aface*T_gas + q_const*Aface;

        end

        % East (outer) radial face between i and i+1
        if i < Nr
            r_face = r_iph;
            Aface = r_face * dtheta * L;
            dE = r_centers(i+1) - ri;
            k_face = harmonic_mean(k_ring(i), k_ring(i+1));
            GE = k_face * Aface / dE;

            idx = idx + 1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + GE;
            q = (j-1)*Nr + (i+1);
            idx = idx + 1; I(idx)=p; J(idx)=q; V(idx)=V(idx) - GE;
        else
            % Outer boundary: convection + radiation to ambient
            Tw0 = Tprev(i,j);
            [h_eff, q_const] = rad_lin_coeff(BC.eps_out, BC.sigma, Tw0, BC.T_amb);
            h_total = BC.h_out + h_eff;

            r_face = r_iph;
            Aface = r_face * dtheta * L;

            idx = idx + 1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + h_total*Aface;
            bP = bP + h_total*Aface*BC.T_amb + q_const*Aface;
        end

        % ---------- Theta conduction neighbors (periodic) ----------
        % FV in theta: Gtheta = k * A_theta / ds, with A_theta = dr_i*L and ds = ri*dtheta
        dr_i = (r_iph - r_imh);
        Atheta = dr_i * L;
        ds = ri * dtheta;

        % Use ring k; if you later want theta-varying k, replace accordingly
        ktheta = k_ring(i);
        Gth = ktheta * Atheta / ds;

        jplus = jp(j);
        jminus = jm(j);

        p_plus  = (jplus-1)*Nr + i;
        p_minus = (jminus-1)*Nr + i;

        idx = idx + 1; I(idx)=p; J(idx)=p;      V(idx)=V(idx) + Gth;
        idx = idx + 1; I(idx)=p; J(idx)=p_plus; V(idx)=V(idx) - Gth;

        idx = idx + 1; I(idx)=p; J(idx)=p;       V(idx)=V(idx) + Gth;
        idx = idx + 1; I(idx)=p; J(idx)=p_minus; V(idx)=V(idx) - Gth;

        % ---------- Transient term ----------
        idx = idx + 1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + aP;

        b(p) = bP;
    end
end

% Trim and build sparse
I = I(1:idx); J = J(1:idx); V = V(1:idx);
A = sparse(I, J, V, N, N);

end

function km = harmonic_mean(k1, k2)
km = 2.0*k1*k2 / (k1 + k2 + eps);
end

function [h_rad, q_const] = rad_lin_coeff(epsr, sigma, Tw_prev, Tenv)
% Linearize eps*sigma*(Tenv^4 - Tw^4) about Tw_prev:
% eps*sigma*(Tenv^4 - Tw^4) ≈ h_rad*(Tenv - Tw) + q_const
% with h_rad = 4*eps*sigma*Tw_prev^3
% and q_const chosen so expression is exact at Tw = Tw_prev.
Tw = max(Tw_prev, 1e-6);
h_rad = 4*epsr*sigma*Tw^3;
q_exact_at_prev = epsr*sigma*(Tenv^4 - Tw^4);
q_const = q_exact_at_prev - h_rad*(Tenv - Tw);
end

function [Tavg, Tmax] = outer_wall_stats(T)
Tout = T(end,:);
Tavg = mean(Tout);
Tmax = max(Tout);
end

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

function E = solid_internal_energy(T, r_faces, dtheta, L, rhoCp_ring)
% Total internal energy proxy: sum(rho*cp*T*V) [J], per unit length L
% T is Nr x Ntheta, rhoCp_ring is Nr x 1
Nr = size(T,1);
Ntheta = size(T,2);

E = 0.0;
for i = 1:Nr
    r_imh = r_faces(i);
    r_iph = r_faces(i+1);
    Vring = 0.5*(r_iph^2 - r_imh^2) * dtheta * L;   % volume per cell in theta
    % sum over theta for ring i
    E = E + rhoCp_ring(i) * Vring * sum(T(i,:));
end
end

function [dQ_in, dQ_out] = boundary_heat_step(Tprev, r_faces, dtheta, L, BC, isBurning, T_gas, dt)
% Compute boundary heat added to the solid over the timestep using the same
% Robin + radiation model as the solver (linearized around Tprev wall).
%
% Returns dQ_in (positive into solid) and dQ_out (positive leaving solid).

Nr = size(Tprev,1);
Ntheta = size(Tprev,2);

% Inner boundary face area per theta cell: A = r_in * dtheta * L
r_in = r_faces(1);
Ain_cell = r_in * dtheta * L;

% Outer boundary face area per theta cell: A = r_out * dtheta * L
r_out = r_faces(end);
Aout_cell = r_out * dtheta * L;

dQ_in = 0.0;
dQ_out = 0.0;

% Inner boundary heat during and after burn: Robin active (per your new model)
% If you still have an adiabatic option, gate it here, but you removed it.
for j = 1:Ntheta
    Tw = Tprev(1,j);
    [h_rad, q_const] = rad_lin_coeff(BC.eps_in, BC.sigma, Tw, T_gas);
    h_tot = BC.h_in + h_rad;

    % q'' into solid (W/m^2): h_tot*(Tgas - Tw) + q_const
    qpp_in = h_tot*(T_gas - Tw) + q_const;
    dQ_in = dQ_in + qpp_in * Ain_cell * dt;
end

% Outer boundary: Robin to ambient, heat leaves solid when Tw > Tamb
for j = 1:Ntheta
    Tw = Tprev(Nr,j);
    [h_rad, q_const] = rad_lin_coeff(BC.eps_out, BC.sigma, Tw, BC.T_amb);
    h_tot = BC.h_out + h_rad;

    % q'' leaving solid (W/m^2): h_tot*(Tw - Tamb) - q_const
    % Because solver uses +h_tot*(Tamb - Tw) + q_const as flux INTO solid.
    % So outward flux = -(into flux) = h_tot*(Tw - Tamb) - q_const.
    qpp_out = h_tot*(Tw - BC.T_amb) - q_const;
    dQ_out = dQ_out + qpp_out * Aout_cell * dt;
end
end

function [dQ_in, dQ_out] = boundary_heat_from_Ab(Tnew, Tprev, r_faces, dtheta, L, BC, T_gas, dt)
% Compute boundary heat using the SAME linearized coefficients as assembly,
% evaluated consistently with the BE solve at Tnew.
%
% dQ_in  = positive heat INTO solid at inner boundary
% dQ_out = positive heat LEAVING solid at outer boundary

Nr = size(Tnew,1);
Ntheta = size(Tnew,2);

r_in  = r_faces(1);
r_out = r_faces(end);

Ain  = r_in  * dtheta * L;
Aout = r_out * dtheta * L;

dQ_in  = 0.0;
dQ_out = 0.0;

% inner boundary (environment = T_gas)
for j = 1:Ntheta
    Tw0 = Tprev(1,j);                % linearization point (matches assembly)
    [h_rad, q_const] = rad_lin_coeff(BC.eps_in, BC.sigma, Tw0, T_gas);
    h_tot = BC.h_in + h_rad;

    K = h_tot * Ain;                 % diagonal contribution
    F = h_tot * Ain * T_gas + q_const * Ain;   % RHS contribution

    q_in = (F - K*Tnew(1,j));         % W (since K is W/K, F is W)
    dQ_in = dQ_in + q_in * dt;        % J
end

% outer boundary (environment = Tamb)
for j = 1:Ntheta
    Tw0 = Tprev(Nr,j);
    [h_rad, q_const] = rad_lin_coeff(BC.eps_out, BC.sigma, Tw0, BC.T_amb);
    h_tot = BC.h_out + h_rad;

    K = h_tot * Aout;
    F = h_tot * Aout * BC.T_amb + q_const * Aout;

    q_into_solid = (F - K*Tnew(Nr,j));    % W into solid (negative when cooling)
    q_out = max(0, -q_into_solid);        % W leaving solid
    dQ_out = dQ_out + q_out * dt;
end
end

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
    figTitle = 'Energy diagnostics (with deletion diagnostics)';
else
    figTitle = 'Energy diagnostics';
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
