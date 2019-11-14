function g = GetRhoU(p,g)

g.rho = sum(g.fIn); 

% original code:
% g.ux = reshape ( (g.cx * reshape(g.fIn,p.Q,p.lx*p.ly)), 1,p.lx,p.ly ) ./g.rho; 
% g.uy = reshape ( (g.cy * reshape(g.fIn,p.Q,p.lx*p.ly)), 1,p.lx,p.ly ) ./g.rho; 

% modified code:
g.ux = reshape ( (g.cx * reshape(g.fIn,p.Q,p.lx*p.ly)), 1,p.lx,p.ly );
g.uy = reshape ( (g.cy * reshape(g.fIn,p.Q,p.lx*p.ly)), 1,p.lx,p.ly ); 

g.rho = squeeze(g.rho);
g.ux = squeeze(g.ux);
g.uy = squeeze(g.uy);
