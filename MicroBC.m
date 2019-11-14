function g = MicroBC(p,g)

% MICROSCOPIC BOUNDARY CONDITIONS 


% %% Surface: constant velocity
% 
% xsurf = 2:(p.lx-1); % x indices to apply surface BC
% ysurf = p.ly; % y index of surface
% 
% % Zou-He BC
% g.fIn(5,xsurf,ysurf) = g.fIn(3,xsurf,ysurf) - 2/3*g.rho(:,xsurf,ysurf).*g.uy(:,xsurf,ysurf); 
% g.fIn(9,xsurf,ysurf) = g.fIn(7,xsurf,ysurf) + 1/2*(g.fIn(4,xsurf,ysurf)-g.fIn(2,xsurf,ysurf))+ ... 
%                 1/2*g.rho(:,xsurf,ysurf).*g.ux(:,xsurf,ysurf) - 1/6*g.rho(:,xsurf,ysurf).*g.uy(:,xsurf,ysurf); 
% g.fIn(8,xsurf,ysurf) = g.fIn(6,xsurf,ysurf) + 1/2*(g.fIn(2,xsurf,ysurf)-g.fIn(4,xsurf,ysurf))- ... 
%                 1/2*g.rho(:,xsurf,ysurf).*g.ux(:,xsurf,ysurf) - 1/6*g.rho(:,xsurf,ysurf).*g.uy(:,xsurf,ysurf);
            

%% Surface: FREE-SLIP (specular reflection)

% Note that we have to use the appropriate reflection indices depending on
% which boundary is free-slip.

% Assign the reflected distributions on fluid nodes to the appropriate
% outgoing directions of solid nodes (this method didn't work)
% g.fOut(sub2ind3D(siz,g.opp(g.fs.dir)',g.fs.xs,g.fs.ys)) = g.fOut(sub2ind3D(siz,g.fs.dist,g.fs.xf,g.fs.yf));


if g.fsW
    g.fOut(g.fsWind,g.fsWrecip) = g.fOut(:,g.fsWdonor);
end

if g.fsE
    g.fOut(g.fsEind,g.fsErecip) = g.fOut(:,g.fsEdonor);
end

if g.fsS
    g.fOut(g.fsSind,g.fsSrecip) = g.fOut(:,g.fsSdonor);
end

if g.fsN
    g.fOut(g.fsNind,g.fsNrecip) = g.fOut(:,g.fsNdonor);
end


% %% Inlet (left boundary)
% 
% in = 1; % x index of inlet
% flin = 2:(p.ly-1);
% % flin = ceil(g.z(in)):p.ly;
% 
% % Zou/He BC
% g.fIn(2,in,flin) = g.fIn(4,in,flin) + 2/3*g.rho(:,in,flin).*g.ux(:,in,flin); 
% g.fIn(6,in,flin) = g.fIn(8,in,flin) + 1/2*(g.fIn(5,in,flin)-g.fIn(3,in,flin)) ... 
%                                 + 1/2*g.rho(:,in,flin).*g.uy(:,in,flin) ...
%                                 + 1/6*g.rho(:,in,flin).*g.ux(:,in,flin); 
% g.fIn(9,in,flin) = g.fIn(7,in,flin) + 1/2*(g.fIn(3,in,flin)-g.fIn(5,in,flin)) ... 
%                                 - 1/2*g.rho(:,in,flin).*g.uy(:,in,flin) ...
%                                 + 1/6*g.rho(:,in,flin).*g.ux(:,in,flin); 

%try boundary condition of zhang kwok 2006



% %% Outlet (right boundary)
% 
% out = p.lx; % x index of inlet
% flout = 2:(p.ly-1);
% % flout = ceil(g.z(out)):p.ly;
% 
% 
% % Zou/He BC
% g.fIn(4,out,flout) = g.fIn(2,out,flout) - 2/3*g.rho(:,out,flout).*g.ux(:,out,flout); 
% g.fIn(8,out,flout) = g.fIn(6,out,flout) + 1/2*(g.fIn(3,out,flout)-g.fIn(5,out,flout)) ... 
%                                   - 1/2*g.rho(:,out,flout).*g.uy(:,out,flout) ...
%                                   - 1/6*g.rho(:,out,flout).*g.ux(:,out,flout); 
% g.fIn(7,out,flout) = g.fIn(9,out,flout) + 1/2*(g.fIn(5,out,flout)-g.fIn(3,out,flout)) ... 
%                                   + 1/2*g.rho(:,out,flout).*g.uy(:,out,flout) ...
%                                   - 1/6*g.rho(:,out,flout).*g.ux(:,out,flout); 

                              

%% Bed: MODIFIED NO-SLIP (Bouzidi)

if any(strcmpi(struct2cell(p.bdy),'bouzidi')),
% the struct g.ns is updated in Bouzidi.m (called by UpdateBC.m) each time the boundary changes

% standard bounce-back for the solid nodes so their velocity remains zero
g.fOut(:,g.bbRegion) = g.fIn(g.opp,g.bbRegion); 

% LHS quantities in Bouzidi et al.'s equation 5, using post-collision and post-streaming f's
fqlthalf = 2*g.ns.q.*g.fOut(g.ns.fd) + (1-2*g.ns.q).*g.fOut(g.ns.f2d);
fqgthalf = (2*g.ns.q-1)./(2*g.ns.q).*g.fOut(g.ns.fu) + 1./(2*g.ns.q).*g.fOut(g.ns.fd);

f(g.ns.q < 0.5) = fqlthalf(g.ns.q < 0.5);
f(g.ns.q >= 0.5) = fqgthalf(g.ns.q >= 0.5);


% Assign the f's to the appropriate outgoing directions of solid nodes
g.fOut(g.ns.fout) = f;
end

if any(strcmpi(struct2cell(p.bdy),'noslip')),
% %% Bed: NO-SLIP ("full-way" bounce-back):
% 
g.fOut(:,g.bbRegion) = g.fIn(g.opp,g.bbRegion); 
end
