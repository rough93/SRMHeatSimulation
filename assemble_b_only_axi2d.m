function b = assemble_b_only_axi2d(Tprev, grid, matField, BC, burn, BCinner, dt, cfg)
%ASSEMBLE_B_ONLY_AXI2D Build RHS vector b for 2D axisymmetric BE solve on solid nodes only.
%
% MUST match assemble_BE_2d_axisym.m exactly in:
%   - unknown ordering (solidList = find(burn.isSolid))
%   - face radii and truncated control volume at the moving interface
%   - which faces get Robin terms:
%       inner Robin when west face borders cavity (or i==1)
%       outer Robin at outer boundary (or if east neighbor is cavity fallback)

Nr = numel(grid.r);
Nz = numel(grid.z);

dr = grid.dr;
dz = grid.dz;

isSolid = burn.isSolid;
solidList = find(isSolid);
Ns = numel(solidList);

b = zeros(Ns, 1);

for idx = 1:Ns
    lin = solidList(idx);
    [i, j] = ind2sub([Nr, Nz], lin);

    ri = grid.r(i);

    % ---- Face radii ----
    rE = ri + dr/2;
    rW = max(ri - dr/2, 0);

    % If west face borders cavity, the west face is the burn front radius
    bordersCavityOnWest = ( (i == 1) || ((i > 1) && ~isSolid(i-1, j)) );
    if bordersCavityOnWest
        rW = max(burn.r_front, 0);
    end

    % ---- FV geometry (2*pi cancels) ----
    V      = 0.5 * (rE^2 - rW^2) * dz;
    Ar_in  = rW * dz;
    Ar_out = rE * dz;

    % ---- Transient RHS ----
    ap  = matField.rhocp(i,j) * V / dt;
    rhs = ap * Tprev(i,j);

    % ---- Inner Robin RHS (matches A-builder) ----
    if bordersCavityOnWest && (BCinner.h > 0)
        rhs = rhs + (BCinner.h * Ar_in) * BCinner.Tinf;
    end

    % ---- Outer Robin RHS (matches A-builder) ----
    bordersCavityOnEast = ( (i < Nr) && ~isSolid(i+1, j) );
    isOuterBoundary     = (i == Nr);

    if (isOuterBoundary || bordersCavityOnEast) && (BC.h_eff > 0)
        rhs = rhs + (BC.h_eff * Ar_out) * BC.T_amb;
    end

    % Axial ends adiabatic => no RHS contribution

    b(idx) = rhs;
end

end
