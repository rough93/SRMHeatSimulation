function [A, b] = assemble_BE_2d_axisym(Tprev, grid, matField, BC, burn, BCinner, dt, cfg)
%ASSEMBLE_BE_2D_AXISYM Backward-Euler FV axisymmetric r-z conduction on SOLID nodes only.
%
% Unknown ordering: solidList = find(burn.isSolid) (column-major), idx=1..Ns
% Governing: rhocp dT/dt = 1/r d/dr (k r dT/dr) + d/dz (k dT/dz)
%
% Inner BC: Robin to gas/cavity applied on cells adjacent to cavity (moving interface)
% Outer BC: Robin to ambient at i==Nr
% Axial ends: adiabatic by default

Nr = numel(grid.r);
Nz = numel(grid.z);

dr = grid.dr;
dz = grid.dz;

isSolid = burn.isSolid;
solidList = find(isSolid);
Ns = numel(solidList);

id = zeros(Nr, Nz);
id(solidList) = 1:Ns;

A = spalloc(Ns, Ns, 7*Ns);
b = zeros(Ns, 1);

harmonic = @(a,b_) (2*a*b_) / (a + b_ + eps);

axial_adiabatic = true;
if isfield(cfg,'axial_adiabatic')
    axial_adiabatic = cfg.axial_adiabatic;
end

for idx = 1:Ns
    lin = solidList(idx);
    [i, j] = ind2sub([Nr, Nz], lin);

    ri = grid.r(i);

    % ---------------- Face radii (axisymmetric FV, 2*pi cancels) ----------------
    rE = ri + dr/2;

    % Default west face at geometric face location
    rW = max(ri - dr/2, 0);

    % If this solid cell borders the cavity on its west face, the west face is the burn front.
    % This is the key fix: truncated control volume at the moving interface.
    bordersCavityOnWest = ( (i == 1) || ((i > 1) && ~isSolid(i-1, j)) );
    if bordersCavityOnWest
        rW = max(burn.r_front, 0);
    end

    % FV geometry (exact, not ri*dr approximations)
    % Volume factor (per unit theta): integral(r dr) * dz = 0.5*(rE^2 - rW^2)*dz
    V  = 0.5 * (rE^2 - rW^2) * dz;

    % Radial face area factors (per unit theta): r_face * dz
    Ar_in  = rW * dz;
    Ar_out = rE * dz;

    % Axial face area factor (per unit theta): integral(r dr) = 0.5*(rE^2 - rW^2)
    Az = 0.5 * (rE^2 - rW^2);

    ap  = matField.rhocp(i,j) * V / dt;
    rhs = ap * Tprev(i,j);

    kC = matField.k(i,j);

    % ---------------- r- (west / inner) ----------------
    if i > 1 && isSolid(i-1, j)
        % Conduction to solid neighbor
        kN = matField.k(i-1, j);
        kf = harmonic(kC, kN);
        a  = kf * Ar_in / dr;

        ap = ap + a;
        A(idx, id(i-1, j)) = A(idx, id(i-1, j)) - a;

    else
        % Cavity interface or i==1 boundary: Robin to gas/cavity
        if BCinner.h > 0
            a   = BCinner.h * Ar_in;
            ap  = ap + a;
            rhs = rhs + a * BCinner.Tinf;
        end
    end

    % ---------------- r+ (east / outer) ----------------
    if i < Nr && isSolid(i+1, j)
        kN = matField.k(i+1, j);
        kf = harmonic(kC, kN);
        a  = kf * Ar_out / dr;

        ap = ap + a;
        A(idx, id(i+1, j)) = A(idx, id(i+1, j)) - a;

    else
        % Outer boundary (or fallback): Robin to ambient
        if BC.h_eff > 0
            a   = BC.h_eff * Ar_out;
            ap  = ap + a;
            rhs = rhs + a * BC.T_amb;
        end
    end

    % ---------------- z- ----------------
    if j > 1 && isSolid(i, j-1)
        kN = matField.k(i, j-1);
        kf = harmonic(kC, kN);
        a  = kf * Az / dz;

        ap = ap + a;
        A(idx, id(i, j-1)) = A(idx, id(i, j-1)) - a;
    else
        % adiabatic unless you implement axial Robin
        if ~axial_adiabatic
            % not implemented
        end
    end

    % ---------------- z+ ----------------
    if j < Nz && isSolid(i, j+1)
        kN = matField.k(i, j+1);
        kf = harmonic(kC, kN);
        a  = kf * Az / dz;

        ap = ap + a;
        A(idx, id(i, j+1)) = A(idx, id(i, j+1)) - a;
    else
        if ~axial_adiabatic
            % not implemented
        end
    end

    A(idx, idx) = A(idx, idx) + ap;
    b(idx) = rhs;
end

end
