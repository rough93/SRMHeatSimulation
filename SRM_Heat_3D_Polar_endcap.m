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

time.dt = 0.050;
time.t_final = 100.0;

cfg.store_every = 10;     % store every N steps
cfg.max_cells = 2e6;      % safety
cfg.Ntheta = 2;          % circumferential resolution (>= 8)
cfg.Nz = 500;              % axial resolution (>=1)
cfg.endBC = 'adiabatic';    % 'adiabatic' or 'robin' at z=0 and z=L

% ---- Optional endcap / bulkhead thermal mass (lumped nodes coupled to z-ends) ----
cfg.endcap.enable = true;        % set false to disable endcap nodes
cfg.endcap.ends   = 'both';      % 'z0', 'zL', or 'both'
cfg.endcap.thickness = 0.012;    % [m] effective bulkhead thickness participating thermally
cfg.endcap.k      = 205;         % [W/m-K] bulkhead conductivity (e.g., aluminum ~205)
cfg.endcap.rho    = 2700;        % [kg/m^3]
cfg.endcap.cp     = 900;         % [J/kg-K]
cfg.endcap.h_in   = 200;         % [W/m^2-K] bulkhead side exposed to hot gas (inside motor)
cfg.endcap.eps_in = 0.8;         % emissivity for inside radiation
cfg.endcap.h_out  = 30;          % [W/m^2-K] bulkhead side exposed to ambient (outside)
cfg.endcap.eps_out = 0.8;        % emissivity for outside radiation

BC.h_z0 = 30;               % end convection at z=0 if endBC='robin'
BC.h_zL = 30;               % end convection at z=L if endBC='robin'
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

Nz = cfg.Nz;
dz = cfg.L / Nz;
%z_centers = ((1:Nz) - 0.5) * dz;

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
E_solid_hist(1) = solid_internal_energy(T, r_faces, dtheta, dz, rhoCp_ring);
%times = [];
%Tout_avg = [];
%Tout_max = [];
snap = struct('t',{},'r_centers',{},'theta',{},'T',{});

for step = 0:nsteps
    t = step * time.dt;

    % ---- store ----
    if mod(step, cfg.store_every)==0
        %times(end+1,1) = t;
        T2 = mean(T,3);  % axial average for plotting/stats
        %[Tavg, Tmax] = outer_wall_stats(T2);
        %Tout_avg(end+1,1) = Tavg;
        %Tout_max(end+1,1) = Tmax;

        snap(end+1).t = t;
        snap(end).r_centers = r_centers;
        snap(end).theta = theta_centers;
        snap(end).T = T2;

        % ---- Store 3D outer-surface temperature map for animation ----
        % Tout3D_store will be Ntheta x Nz x Nt_store
        if ~exist('Tout3D_store','var')
            Tout3D_store = [];
            t_store = [];
        end

        t_store(end+1,1) = t;

        % outer surface temperature as a function of theta and z
        Tout_theta_z = squeeze(T(end,:,:));        % size: [Ntheta x Nz]
        Tout3D_store(:,:,end+1) = Tout_theta_z;    % append along 3rd dim


        % 1D equivalents for original-code style plots
        snap(end).Tbar_r = mean(T2, 2);
        snap(end).Tmax_r = max(T2, [], 2);
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
    [A, b] = assemble_BE_polar_3D( ...
        Tprev, Tcap0, TcapL, r_faces, r_centers, dtheta, dz, ...
        k_ring, rhoCp_ring, time.dt, ...
        BC, isBurning, T_gas, cfg);
    asm_time = toc(t_asm);

    t_solve = tic;
    Tvec = A \ b;
    solve_time = toc(t_solve);

    Nwall = Nr * Ntheta * Nz;

    % Wall temperatures live in the first Nwall entries
    T = reshape(Tvec(1:Nwall), [Nr, Ntheta, Nz]);

    % Extract endcap temperatures (if present)
    if isfield(cfg,'endcap') && isfield(cfg.endcap,'enable') && cfg.endcap.enable
        if ~isfield(cfg.endcap,'ends'); cfg.endcap.ends = 'both'; end

        if strcmpi(cfg.endcap.ends,'both')
            Tcap0 = Tvec(Nwall + 1);
            TcapL = Tvec(Nwall + 2);

        elseif strcmpi(cfg.endcap.ends,'z0')
            Tcap0 = Tvec(Nwall + 1);
            TcapL = [];

        elseif strcmpi(cfg.endcap.ends,'zL')
            Tcap0 = [];
            TcapL = Tvec(Nwall + 1);
        end
    end


    % Extract endcap temperatures (if present)
    if isfield(cfg,'endcap') && isfield(cfg.endcap,'enable') && cfg.endcap.enable
        Nwall = Nr*Ntheta*Nz;
        if strcmpi(cfg.endcap.ends,'both')
            Tcap0 = Tvec(Nwall+1);
            TcapL = Tvec(Nwall+2);
        elseif strcmpi(cfg.endcap.ends,'z0')
            Tcap0 = Tvec(Nwall+1);
            TcapL = [];
        elseif strcmpi(cfg.endcap.ends,'zL')
            Tcap0 = [];
            TcapL = Tvec(Nwall+1);
        end
    end
    step_time = toc(t_start_step);

    % =========================
    % Boundary heat (consistent with assembler: uses Tprev for rad linearization)
    % =========================
    dQ_in  = 0.0;   % positive INTO solid
    dQ_out = 0.0;   % positive LEAVING solid

    % Inner/outer areas per cell use dz (not cfg.L), because we explicitly resolve Nz slices.
    for kk = 1:Nz
        for j = 1:Ntheta
            % ---------- INNER boundary (i=1) ----------
            Tw0 = Tprev(1,j,kk);
            [h_eff, q_const] = rad_lin_coeff(BC.eps_in, BC.sigma, Tw0, T_gas);
            h_total = BC.h_in + h_eff;

            Aface_in = r_faces(1) * dtheta * dz;   % area of inner radial face for this (theta,z) cell
            dQ_in = dQ_in + time.dt * ( h_total*Aface_in*(T_gas - T(1,j,kk)) + q_const*Aface_in );

            % ---------- OUTER boundary (i=Nr) ----------
            Tw0 = Tprev(Nr,j,kk);
            [h_eff, q_const] = rad_lin_coeff(BC.eps_out, BC.sigma, Tw0, BC.T_amb);
            h_total = BC.h_out + h_eff;

            Aface_out = r_faces(end) * dtheta * dz; % outer radial face at r_{Nr+1/2}
            dQ_into_solid_outer = time.dt * ( h_total*Aface_out*(BC.T_amb - T(Nr,j,kk)) + q_const*Aface_out );

            % bookkeeping: define dQ_out as positive LEAVING solid
            dQ_out = dQ_out - dQ_into_solid_outer;
        end
    end

    dQnet = dQ_in - dQ_out;


    % ---- Energy balance for this step (BEFORE deleting rings) ----
    % --- Energies computed on the same domain as the solve (pre-deletion) ---
    E_prev = E_solid_hist(step+1);  % energy of the carried state at start of this step
    E_new  = solid_internal_energy(T, r_faces, dtheta, dz, rhoCp_ring);

    % Energy that leaves the computational domain if we delete ring 1 AFTER the solve
    E_removed = 0.0;
    if willDelete
        r_imh = r_faces(1);
        r_iph = r_faces(2);
        Vcell = 0.5*(r_iph^2 - r_imh^2) * dtheta * dz;   % per-(theta,z)-cell volume
        E_removed = rhoCp_ring(1) * Vcell * sum(T(1,:,:),'all');    % sum over theta
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

    % Only meaningful if you are about to delete ring 1
    if willDelete && Nr >= 2
        % We want the conduction flow between ring 2 and ring 1:
        % Q = dt * sum_over_theta,z( G12 * (T2 - T1) )
        %
        % For BE matrix, the coupling entry A(p1,p2) is -G12 (for that cell),
        % so G12 = -A(p1,p2).

        % Build indices for all (j,kk) pairs
        [Jgrid, Kgrid] = ndgrid(1:Ntheta, 1:Nz);
        Jlist = Jgrid(:);
        Klist = Kgrid(:);

        p1 = sub2ind([Nr, Ntheta, Nz], ones(size(Jlist)), Jlist, Klist);   % (i=1, j, k)
        p2 = sub2ind([Nr, Ntheta, Nz], 2*ones(size(Jlist)), Jlist, Klist); % (i=2, j, k)

        % Extract the coupling coefficients A(p1,p2)
        lin = sub2ind(size(A), p1, p2);
        G12 = full(-A(lin));   % length = Ntheta*Nz

        % Temperature difference (T2 - T1) for those same cells
        dT12 = reshape(T(2,:,:) - T(1,:,:), [], 1);  % length = Ntheta*Nz

        Qcond12_fromA = time.dt * sum(G12 .* dT12);
    end

    Qcond12A_hist(step+1) = Qcond12_fromA;


    % ---- status output ----
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
        T(1,:,:) = [];
        r_centers(1) = [];
        r_faces(1) = [];
        layer_idx_of_ring(1) = [];

        Nr = numel(r_centers);
        [k_ring, rhoCp_ring] = ring_props(layers, layer_idx_of_ring);
    end

end

plot_results(snap, layers, geom, burn, BC);

%% ==================== 3D tube animation (outer surface) ====================
% Requires: Tout3D_store (Ntheta x Nz x Nt_store), t_store, theta_centers, z_centers, r_faces

%% ==================== 3D tube animation + GIF export (outer surface) ====================
% Requires: Tout3D_store (Ntheta x Nz x Nt_store), t_store, theta_centers, r_faces

if exist('Tout3D_store','var') && ~isempty(Tout3D_store)
    r_outer = r_faces(end);
    z_centers = ((1:cfg.Nz) - 0.5) * (cfg.L / cfg.Nz);

    % Build cylinder surface grid
    [TH, ZZ] = meshgrid(theta_centers, z_centers);
    XX = r_outer * cos(TH);
    YY = r_outer * sin(TH);

    % Initial color data
    C0 = Tout3D_store(:,:,1).';   % [Nz x Ntheta]

    fig = figure('Name','3D Outer Surface Temperature', 'Color','w');
    h = surf(XX, YY, ZZ, C0, 'EdgeColor','none');
    axis equal tight;
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    cb = colorbar; ylabel(cb, 'Temperature [K]');
    colormap(turbo);
    view(35, 20);
    camlight headlight; lighting gouraud;

    % Fix color scaling for entire animation
    Tmin = min(Tout3D_store, [], 'all');
    Tmax = max(Tout3D_store, [], 'all');
    caxis([Tmin Tmax]);

    % Frame limits (protect against mismatch)
    NtT = size(Tout3D_store, 3);
    Ntt = numel(t_store);
    Nt  = min(NtT, Ntt);

    % GIF settings
    gifname = 'outer_wall_temperature.gif';
    gif_dt  = 0.05;   % seconds per frame (20 FPS)

    % --------- Animation + GIF write ---------
    for kk = 1:Nt
        Ck = Tout3D_store(:,:,kk).';
        set(h, 'CData', Ck);
        title(sprintf('Outer wall T, t = %.2f s', t_store(kk)));
        drawnow;

        % Capture frame
        frame = getframe(fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        if kk == 1
            imwrite(imind, cm, gifname, 'gif', ...
                'Loopcount', inf, 'DelayTime', gif_dt);
        else
            imwrite(imind, cm, gifname, 'gif', ...
                'WriteMode', 'append', 'DelayTime', gif_dt);
        end
    end

    fprintf('GIF written to %s (%d frames)\n', gifname, Nt);

else
    warning('Tout3D_store not found or empty. Add the store block inside the sim loop first.');
end




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
    k_ring, rhoCp_ring, dt, BC, ~, T_gas, cfg)
% Assemble sparse Backward Euler FV system for 3D polar conduction (r,theta,z) with optional lumped endcaps.
% Unknown ordering for wall cells: p = sub2ind([Nr, Ntheta, Nz], i, j, k)
% Optional endcap nodes (lumped):
%   If enabled, extra unknown(s) appended after wall cells: [Tcap0; TcapL] depending on cfg.endcap.ends.

Nr = numel(r_centers);
Ntheta = size(Tprev,2);
Nz = size(Tprev,3);
Nwall = Nr * Ntheta * Nz;

useEndcap = isfield(cfg,'endcap') && isfield(cfg.endcap,'enable') && cfg.endcap.enable;
cap_z0 = false; cap_zL = false;
if useEndcap
    if ~isfield(cfg.endcap,'ends'); cfg.endcap.ends = 'both'; end
    cap_z0 = strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'z0');
    cap_zL = strcmpi(cfg.endcap.ends,'both') || strcmpi(cfg.endcap.ends,'zL');
end
Ncap = double(cap_z0) + double(cap_zL);

N = Nwall + Ncap;

% Preallocate sparse
nz_est = N * 14;
I = zeros(nz_est,1);
J = zeros(nz_est,1);
V = zeros(nz_est,1);
b = zeros(N,1);
idx = 0;

% Convenience: periodic theta indexing
jp = @(j) (j < Ntheta) * (j+1) + (j==Ntheta) * 1;
jm = @(j) (j > 1)     * (j-1) + (j==1)     * Ntheta;

% Endcap node indices (if present)
p_cap0 = [];
p_capL = [];
if cap_z0
    p_cap0 = Nwall + 1;
end
if cap_zL
    p_capL = Nwall + 1 + double(cap_z0);
end

% Annulus area for endcap coupling and convection/radiation (matches the wall domain)
r_inner = r_faces(1);
r_outer = r_faces(end);
A_annulus = pi * (r_outer^2 - r_inner^2);     % [m^2]

% Endcap thermal mass
if useEndcap
    if ~isfield(cfg.endcap,'thickness'); cfg.endcap.thickness = dz; end
    if ~isfield(cfg.endcap,'k');  cfg.endcap.k  = 1.0; end
    if ~isfield(cfg.endcap,'rho'); cfg.endcap.rho = 1000; end
    if ~isfield(cfg.endcap,'cp');  cfg.endcap.cp  = 1000; end
    rhoCp_cap = cfg.endcap.rho * cfg.endcap.cp;

    Vcap = A_annulus * cfg.endcap.thickness;  % [m^3]
    aPcap = rhoCp_cap * Vcap / dt;            % [W/K] in BE form

    % Helper to assemble a single endcap node equation

    % Assemble endcap node equations (lumped) if enabled
    if cap_z0
        pcap = p_cap0; Tcap_prev = Tcap0_prev;
        % transient term
        idx = idx + 1; I(idx)=pcap; J(idx)=pcap; V(idx)=V(idx) + aPcap;
        b(pcap) = b(pcap) + aPcap * Tcap_prev;

        % inside (to gas) convection+radiation
        if ~isfield(cfg.endcap,'h_in'); cfg.endcap.h_in = 0; end
        if ~isfield(cfg.endcap,'eps_in'); cfg.endcap.eps_in = 0; end
        Tw0 = Tcap_prev;
        [h_rad_in, q_const_in] = rad_lin_coeff(cfg.endcap.eps_in, BC.sigma, Tw0, T_gas);
        htot_in = cfg.endcap.h_in + h_rad_in;
        if htot_in > 0
            idx = idx + 1; I(idx)=pcap; J(idx)=pcap; V(idx)=V(idx) + htot_in*A_annulus;
            b(pcap) = b(pcap) + htot_in*A_annulus*T_gas + q_const_in*A_annulus;
        end

        % outside (to ambient) convection+radiation
        if ~isfield(cfg.endcap,'h_out'); cfg.endcap.h_out = 0; end
        if ~isfield(cfg.endcap,'eps_out'); cfg.endcap.eps_out = 0; end
        Tw1 = Tcap_prev;
        [h_rad_out, q_const_out] = rad_lin_coeff(cfg.endcap.eps_out, BC.sigma, Tw1, BC.T_amb);
        htot_out = cfg.endcap.h_out + h_rad_out;
        if htot_out > 0
            idx = idx + 1; I(idx)=pcap; J(idx)=pcap; V(idx)=V(idx) + htot_out*A_annulus;
            b(pcap) = b(pcap) + htot_out*A_annulus*BC.T_amb + q_const_out*A_annulus;
        end
    end

    if cap_zL
        pcap = p_capL; Tcap_prev = TcapL_prev;
        % transient term
        idx = idx + 1; I(idx)=pcap; J(idx)=pcap; V(idx)=V(idx) + aPcap;
        b(pcap) = b(pcap) + aPcap * Tcap_prev;

        % inside (to gas) convection+radiation
        if ~isfield(cfg.endcap,'h_in'); cfg.endcap.h_in = 0; end
        if ~isfield(cfg.endcap,'eps_in'); cfg.endcap.eps_in = 0; end
        Tw0 = Tcap_prev;
        [h_rad_in, q_const_in] = rad_lin_coeff(cfg.endcap.eps_in, BC.sigma, Tw0, T_gas);
        htot_in = cfg.endcap.h_in + h_rad_in;
        if htot_in > 0
            idx = idx + 1; I(idx)=pcap; J(idx)=pcap; V(idx)=V(idx) + htot_in*A_annulus;
            b(pcap) = b(pcap) + htot_in*A_annulus*T_gas + q_const_in*A_annulus;
        end

        % outside (to ambient) convection+radiation
        if ~isfield(cfg.endcap,'h_out'); cfg.endcap.h_out = 0; end
        if ~isfield(cfg.endcap,'eps_out'); cfg.endcap.eps_out = 0; end
        Tw1 = Tcap_prev;
        [h_rad_out, q_const_out] = rad_lin_coeff(cfg.endcap.eps_out, BC.sigma, Tw1, BC.T_amb);
        htot_out = cfg.endcap.h_out + h_rad_out;
        if htot_out > 0
            idx = idx + 1; I(idx)=pcap; J(idx)=pcap; V(idx)=V(idx) + htot_out*A_annulus;
            b(pcap) = b(pcap) + htot_out*A_annulus*BC.T_amb + q_const_out*A_annulus;
        end
    end
end

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

            % ---------- Radial conduction neighbors ----------
            if i > 1
                r_face = r_imh;
                Aface  = r_face * dtheta * dz;
                dW     = ri - r_centers(i-1);
                k_face = harmonic_mean(k_ring(i-1), k_ring(i));
                GW     = k_face * Aface / dW;

                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + GW;
                q = sub2ind([Nr,Ntheta,Nz], i-1, j, kk);
                idx=idx+1; I(idx)=p; J(idx)=q; V(idx)=V(idx) - GW;
            else
                % Inner boundary: convection + radiation to gas
                Tw0 = Tprev(i,j,kk);
                [h_eff, q_const] = rad_lin_coeff(BC.eps_in, BC.sigma, Tw0, T_gas);
                h_total = BC.h_in + h_eff;

                r_face = r_imh;
                Aface  = r_face * dtheta * dz;

                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + h_total*Aface;
                bP = bP + h_total*Aface*T_gas + q_const*Aface;
            end

            if i < Nr
                r_face = r_iph;
                Aface  = r_face * dtheta * dz;
                dE     = r_centers(i+1) - ri;
                k_face = harmonic_mean(k_ring(i), k_ring(i+1));
                GE     = k_face * Aface / dE;

                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + GE;
                q = sub2ind([Nr,Ntheta,Nz], i+1, j, kk);
                idx=idx+1; I(idx)=p; J(idx)=q; V(idx)=V(idx) - GE;
            else
                % Outer boundary: convection + radiation to ambient
                Tw0 = Tprev(i,j,kk);
                [h_eff, q_const] = rad_lin_coeff(BC.eps_out, BC.sigma, Tw0, BC.T_amb);
                h_total = BC.h_out + h_eff;

                r_face = r_iph;
                Aface  = r_face * dtheta * dz;

                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + h_total*Aface;
                bP = bP + h_total*Aface*BC.T_amb + q_const*Aface;
            end

            % ---------- Theta conduction neighbors (periodic) ----------
            dr_i = (r_iph - r_imh);
            Atheta = dr_i * dz;
            ds = ri * dtheta;
            Gth = k_ring(i) * Atheta / ds;

            jE = jp(j); jW = jm(j);
            qE = sub2ind([Nr,Ntheta,Nz], i, jE, kk);
            qW = sub2ind([Nr,Ntheta,Nz], i, jW, kk);

            idx=idx+1; I(idx)=p; J(idx)=p;  V(idx)=V(idx) + 2*Gth;
            idx=idx+1; I(idx)=p; J(idx)=qE; V(idx)=V(idx) - Gth;
            idx=idx+1; I(idx)=p; J(idx)=qW; V(idx)=V(idx) - Gth;

            % ---------- Axial conduction neighbors ----------
            Az = 0.5*(r_iph^2 - r_imh^2) * dtheta;
            Gz = k_ring(i) * Az / dz;

            if kk > 1
                q = sub2ind([Nr,Ntheta,Nz], i, j, kk-1);
                idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + Gz;
                idx=idx+1; I(idx)=p; J(idx)=q; V(idx)=V(idx) - Gz;
            else
                % z=0 end: either couple to endcap, or apply endBC
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
                q = sub2ind([Nr,Ntheta,Nz], i, j, kk+1);
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

            % ---------- Transient term ----------
            idx=idx+1; I(idx)=p; J(idx)=p; V(idx)=V(idx) + aP;
            b(p) = bP;
        end
    end
end

% finalize sparse
I = I(1:idx); J = J(1:idx); V = V(1:idx);
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

function [h_rad, q_const] = rad_lin_coeff(eps, sigma, Tw0, Tenv)
%RAD_LIN_COEFF Linearize radiation heat flux about Tw0.
%
% Radiation heat flux into the surface (positive into solid):
%   q'' = eps*sigma*(Tenv^4 - Tw^4)
%
% Linearization about Tw0:
%   q'' ≈ h_rad*(Tenv - Tw) + q_const
% where:
%   h_rad  = 4*eps*sigma*Tw0^3
%   q_const = eps*sigma*(Tenv^4 - Tw0^4) - h_rad*(Tenv - Tw0)

Tw0 = max(Tw0, 1e-6);  % avoid zero/negative issues
h_rad = 4 * eps * sigma * (Tw0^3);
q_const = eps * sigma * (Tenv^4 - Tw0^4) - h_rad * (Tenv - Tw0);
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