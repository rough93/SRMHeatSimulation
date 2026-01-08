function mat = material_map_r(grid, materials, layers, geom)

Nr = numel(grid.r);
Nz = numel(grid.z);

k     = zeros(Nr, Nz);
rhocp = zeros(Nr, Nz);
layerId = zeros(Nr, Nz);

% interface radii (absolute r)
r_if = zeros(numel(layers)-1,1);
acc = geom.r_inner0;
for i = 1:numel(layers)-1
    acc = acc + layers(i).thickness;
    r_if(i) = acc;
end

% absolute boundaries
r_if_prop_liner = r_if(1);
r_if_liner_case = r_if(2);

R = grid.R;

mask_prop  = (R <= r_if_prop_liner);
mask_liner = (R >  r_if_prop_liner) & (R <= r_if_liner_case);
mask_case  = (R >  r_if_liner_case);

% map key -> material
prop  = materials.prop;
liner = materials.liner;
caseM = materials.case;

k(mask_prop)      = prop.k;
rhocp(mask_prop)  = prop.rho * prop.cp;
layerId(mask_prop)= 1;

k(mask_liner)      = liner.k;
rhocp(mask_liner)  = liner.rho * liner.cp;
layerId(mask_liner)= 2;

k(mask_case)      = caseM.k;
rhocp(mask_case)  = caseM.rho * caseM.cp;
layerId(mask_case)= 3;

mat = struct();
mat.k = k;
mat.rhocp = rhocp;
mat.layerId = layerId;

% also store interfaces in r_rel coordinates for plotting
mat.interfaces_rel = [layers(1).thickness, layers(1).thickness + layers(2).thickness];

end
