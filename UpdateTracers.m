function g = UpdateTracers(p,g,ux,uy)


% advect tracers. LBM time step in lattice units is 1. Could use a midpoint
% method for better accuracy, but that would require another interpolation.

dt = p.saveint; % time step in lattice units (number of iterations)

% % Option 1: Forward Euler
% 
% % interpolate velocities at tracer locations
% [uxt,uyt] = InterpTracerVelocity(p,ux,uy,g.xt,g.yt);
% g.xt = g.xt + uxt/p.u0*dt; % note that we convert velocity to lattice units
% g.yt = g.yt + uyt/p.u0*dt;
% 
% % periodic boundaries, even if model didn't have periodic boundaries
% g.xt = mod( g.xt - 0.5, p.lx ) + 0.5;
% g.yt = mod( g.yt - 0.5, p.ly ) + 0.5;

% Option 2: Midpoint method

[uxt,uyt] = InterpTracerVelocity(p,ux,uy,g.xt,g.yt);
xthalf = g.xt + uxt/p.u0*0.5*dt; % note that we convert velocity to lattice units
ythalf = g.yt + uyt/p.u0*0.5*dt;

% periodic boundaries, even if model didn't have periodic boundaries
xthalf = mod( xthalf - 0.5, p.lx ) + 0.5;
ythalf = mod( ythalf - 0.5, p.ly ) + 0.5;

[uxt,uyt] = InterpTracerVelocity(p,ux,uy,xthalf,ythalf); % velocity at half step
g.xt = g.xt + uxt/p.u0*dt; 
g.yt = g.yt + uyt/p.u0*dt;

g.xt = mod( g.xt - 0.5, p.lx ) + 0.5;
g.yt = mod( g.yt - 0.5, p.ly ) + 0.5;


% Eliminate tracers whose closest lattice site is solid

% round coordinates to the nearest lattice site
xtr = round(g.xt);
ytr = round(g.yt);

% convert to linear indices
itr = sub2ind([p.lx p.ly],xtr,ytr);

g.xt = g.xt(~p.obst(itr));
g.yt = g.yt(~p.obst(itr));


%% subfunctions

function [uxs,uys] = InterpTracerVelocity(p,ux,uy,xs,ys,order)

if nargin < 6
    order = 1;
end

% interpolate the velocities uxs,uys at the points xs,zs

% Note that we have already done three things to pre-process the velocity fields:
% 1. force obstacles to have zero velocity 
% 2. force free-slip boundaries to have same velocity as neighboring row/column
% 3. pad velocity arrays according to boundary conditions 

[Y,X] = meshgrid(0:p.ly+1,0:p.lx+1); % note the extension of x and y to allow interpolation at boundaries

% interpolate
if order == 1
    uxs = interp2(Y,X,ux,ys,xs,'*linear');
    uys = interp2(Y,X,uy,ys,xs,'*linear');
else
    uxs = interp2(Y,X,ux,ys,xs,'*cubic');
    uys = interp2(Y,X,uy,ys,xs,'*cubic');
end

