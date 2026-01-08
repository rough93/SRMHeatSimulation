function burn = burn_front_rzt(grid, burnCfg, t)
% Returns burn state for radial-only burn on a fixed r-z grid.

Nr = numel(grid.r);
Nz = numel(grid.z);

burn.isBurning = (t <= burnCfg.t_burn);

burn.r_front = grid.r(1) + burnCfg.rburn * min(t, burnCfg.t_burn);

% first solid radial index
iFront = find(grid.r >= burn.r_front, 1, 'first');
if isempty(iFront)
    iFront = Nr;
end
burn.iFront = iFront;

% solid is everything at and outside the front (radial-only)
burn.isSolid = false(Nr, Nz);
burn.isSolid(iFront:Nr, :) = true;

% inner face nodes are the first solid ring
burn.isInnerFace = false(Nr, Nz);
burn.isInnerFace(iFront, :) = true;

end