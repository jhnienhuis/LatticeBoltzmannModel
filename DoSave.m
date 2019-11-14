function g = DoSave(p,g)
 
% compute and save dimensional versions of requested quantities. Also save
% dimensional time and iteration count.

i = g.n/p.saveint + 1; % This is the ith save. The +1 is because we also saved the initial state

g.iteration(i) = g.n;
g.time(i) = g.n*p.t0;

if any(ismember(p.savequant,'velocity'))
    
    g = GetRhoU(p,g);
    g.velocity.ux(:,:,i) = g.ux*p.u0; % x velocity in m/s
    g.velocity.uy(:,:,i) = g.uy*p.u0; % y velocity in m/s
    
end

if any(ismember(p.savequant,'shearstress'))
    
    if ~any(ismember(p.savequant,'velocity')) % if we aren't saving velocity, we need to calculate it here
        g = GetRhoU(p,g);
    end
    g = GetStress(p,g);
    g.shearstress(:,:,i) = g.tau*p.p0; % shear stress in Pa
    
end

if any(ismember(p.savequant,'pressure'))
    
    if ~any(ismember(p.savequant,'velocity')) % if we aren't saving velocity, we need to calculate it here
        g = GetRhoU(p,g);
    end
    g.pressure(:,:,i) = (1/3)*g.rho*p.p0; % pressure in Pa
        
end

if any(ismember(p.savequant,'bedshearstress'))
    
    if ~any(ismember(p.savequant,'velocity')) % if we aren't saving velocity, we need to calculate it here
        g = GetRhoU(p,g);
    end
    g = GetTauB(p,g);
    g.bedshearstress(:,1,i) = g.TauB*p.p0; % bed shear stress in Pa
    
end