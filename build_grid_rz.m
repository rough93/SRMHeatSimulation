function grid = build_grid_rz(geom, layers, num)
t_wall = sum([layers.thickness]);

r = linspace(geom.r_inner0, geom.r_inner0 + t_wall, num.Nr).';
z = linspace(0, geom.Lz, num.Nz);

dr = r(2) - r(1);
dz = z(2) - z(1);

[R, Z] = ndgrid(r, z);

grid = struct();
grid.r = r;
grid.z = z;
grid.R = R;
grid.Z = Z;
grid.dr = dr;
grid.dz = dz;
grid.t_wall = t_wall;
grid.r_outer = geom.r_inner0 + t_wall;
end
