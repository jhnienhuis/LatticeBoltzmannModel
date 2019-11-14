function g = GetStress(p,g)

% Deviatoric stress tensor
g.s = GetStressTensor(g.omega,p.nu,g.fIn,g.fEq,p.D,g.Scoeff);

% Shear stress
g.ux = reshape(g.ux,p.lx,p.ly);
g.uy = reshape(g.uy,p.lx,p.ly);
g.tau = ShearStress(g.s,g.ux,g.uy);
