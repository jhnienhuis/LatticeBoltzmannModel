function DoPlot(p,g)


switch p.plotquant
    
    case 'velocity'
        
        g = GetRhoU(p,g);
        g.ux = reshape(g.ux,p.lx,p.ly);
        g.uy = reshape(g.uy,p.lx,p.ly);
        u = sqrt(g.ux.*g.ux + g.uy.*g.uy); % velocity magnitude in LBM units
        quantity = u*p.u0; % velocity magnitude in m/s
        
    case 'pressure'
        
        g = GetRhoU(p,g);
        pressure = (1/3)*g.rho; % pressure in LBM units
        quantity = pressure*p.p0; % pressure in Pa
        
    case 'shearstress'
        
        g = GetRhoU(p,g);
        g = GetStress(p,g);
        quantity = g.tau*p.p0; % shear stress in Pa
        
end

DrawLBMPlot(p,g,quantity);
