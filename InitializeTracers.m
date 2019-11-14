function g = InitializeTracers(p)

% initial locations
[g.yt,g.xt] = meshgrid(1:p.tracerSpacing:p.ly,1:p.tracerSpacing:p.lx);

g.xt = g.xt(:);
g.yt = g.yt(:);

g.yt = g.yt + p.tracerSpacing*(rand(size(g.yt))-0.5);
g.xt = g.xt + p.tracerSpacing*(rand(size(g.xt))-0.5);


% periodic boundaries, even if model didn't have periodic boundaries: if we
% lose a tracer off one side, a new one will pop up on the other side,
% conserving the total number of tracers.
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
