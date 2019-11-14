function g = GetTauB(p,g,order)

if nargin < 3
    order = 1;
end

% find angle of surface tangent at each point
theta = atan2(FastCircShift(p.bed,[-1 0]) - FastCircShift(p.bed,[1 0]) , 2);

% unit vector normal to bed
nx = sin(-theta);
ny = cos(-theta);

% unit vector parallel to the bed
px = cos(theta);
py = sin(theta);

% choose distance normal to surface over which we'll evaluate velocity
% gradient by finite difference
dn = 2; % in lattice units

% calculate coordinates of measurement points
xs = (1:p.lx)' + dn*nx;
zs = p.bed + dn*ny;
xs2 = (1:p.lx)' + 2*dn*nx;
zs2 = p.bed + 2*dn*ny;

% interpolate velocity components at measurement points
[uxs uys] = InterpVelocity(p,g,xs,zs,order);
%[uxs2 uys2] = InterpVelocity(p,g,xs2,zs2,order);

u = dot([uxs uys],[px py],2); % % components of interpolated velocity vectors parallel to bed. To transform into vectors for plotting, use quiver(xs,zs,u.*px,u.*py)
% v = dot([uxs uys],[nx ny],2); % components normal to bed, if desired
%u2 = dot([uxs2 uys2],[px py],2);

% boundary stress is rho*nu*du/dz
% Technically, we should interpolate rho as well, but we would need to know it at the boundary, which we do not. 
% Everywhere on the grid, it is within ~0.005 of 1, so it is a reasonable approximation to assume it is 1.
rho = 1;
% if order == 1
      uz = u/dn;
% elseif order == 2
% %    uz = u/dn - (u2-2*u+0)/dn^2*dn;
%     uz = 2*u - 0.5*u2; % assuming dn =1 
% elseif order == 3
%     g = GetStress(p,g);
%     taub = g.tau(  sub2ind( size(g.tau), 1:p.lx, ceil(g.z)' ) );
%     indz1 = sub2ind(size(g.tau),1:p.lx,ceil(g.z'));
%     indz2 = sub2ind(size(g.tau),1:p.lx,ceil(g.z')+1);
%     %taub = g.tau(indz1) - (ceil(g.z')-g.z').*(g.tau(indz2)-g.tau(indz1));
%     %    keyboard;
%     return
% else
%     error('bleh');
% end
g.TauB = rho*p.nu .*uz; % make sure we get the sign right! make u inherit sign of ux above?


% % smooth the taub vector to eliminate pixel-scale noise
%w = p.smoothtaub; % must be odd
%if w
%   taub = smooth( padarray(taub,[(w-1)/2 0],'circular') , w);
%  % taub = medfilt1(padarray(taub,[(w-1)/2 0],'circular'),w);
%  taub = taub((w-1)/2 + (1:p.lx));
end
