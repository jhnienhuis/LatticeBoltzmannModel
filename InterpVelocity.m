function [uxs uys] = InterpVelocity(p,g,xs,ys,order)

if nargin < 5
    order = 1;
end
% interpolate the velocities uxs,uys atthe points xs,zs

uxs = zeros(size(xs));
uys = zeros(size(ys));

% pad periodically in x dir
% ux = padarray(reshape(g.ux,p.lx,p.ly),[1 0],'circular');
% uy = padarray(reshape(g.uy,p.lx,p.ly),[1 0],'circular');
ux = reshape(g.ux,p.lx,p.ly);
uy = reshape(g.uy,p.lx,p.ly);
ux(g.bbRegion) = 0;
uy(g.bbRegion) = 0;
% ux(g.bbRegion) = 0;
% uy(g.bbRegion) = 0;
ux = FastPadPeriodic(ux,[1 0]);
uy = FastPadPeriodic(uy,[1 0]);
if isfield(p,'smoothu') && length(p.smoothu) > 0
    ux = imfilter(ux,p.smoothu,'circular');
    uy = imfilter(uy,p.smoothu,'circular');
end
[Y X] = meshgrid(1:p.ly,0:p.lx+1); % note the extension of x to allow interpolation at periodic boundaries

% separate points surrounded by regular lattice sites, for which we can use
% interp2 (table lookup), from points near the boundary, where we must use
% a method for irregularly spaced data.
% NOTE THAT WE ASSUME NO OVERHANGING BOUNDARIES IN THIS VERSION (i.e., z(x) is singly valued)

% matrix coordinates (adding 1 to the xs because of the padding) of
% bounding lattice sites
xup = ceil(xs)+1;
xdown = floor(xs)+1;
yup = ceil(ys);
ydown = floor(ys);

% linear indices of neighbor nodes
siz = [p.lx+2 p.ly];
ul = sub2ind2D(siz,xdown,yup);
ur = sub2ind2D(siz,xup,yup);
ll = sub2ind2D(siz,xdown,ydown);
lr = sub2ind2D(siz,xup,ydown);


% logical matrix of solid nodes
solid = false(p.lx,p.ly);
solid(g.bbRegion) = true;
% solid = padarray(solid,[1 0],'circular');
solid = FastPadPeriodic(solid,[1 0]);


% points that have LL and/or LR neighbor below the boundary
irr = solid(ll) | solid(lr);

% make lists of points that do & don't require irregular interpolation
xirr = xs(irr);
yirr = ys(irr);
xreg = xs(~irr);
yreg = ys(~irr);


% for the regular points, use interp2
if ~isempty(xreg)
    if order == 1
        uxs(~irr) = interp2(Y,X,ux,yreg,xreg,'*linear');
        uys(~irr) = interp2(Y,X,uy,yreg,xreg,'*linear');
        % uxs(~irr) = interp2(Y,X,ux,yreg,xreg);
        % uys(~irr) = interp2(Y,X,uy,yreg,xreg);
    else
        uxs(~irr) = interp2(Y,X,ux,yreg,xreg,'*cubic');
        uys(~irr) = interp2(Y,X,uy,yreg,xreg,'*cubic');
    end
end


% for the irregular points, use TriScatteredInterp (Delaunay Triangulation with linear interpolation) 
if ~isempty(xirr)
    
    % make lists of the known points. we include their fluid neighbors and 
    % the boundary points at the neighboring x nodes
    uselattice = false(size(solid));
    uselattice([ul(irr); ur(irr); ll(irr); lr(irr)]) = true; % eliminate points that are not near boundary
    uselattice(solid) = false; % eliminate neighboring points that are not fluid
    
    xb = (0:p.lx+1)';
%     zb = padarray(p.bed,[1 0],'circular');
    zb = FastPadPeriodic(p.bed,[1 0]);
    usebdy = false(size(xb));
    usebdy([xdown(irr); xup(irr)]) = true;
    uxb = zeros(size(xb));
    uyb = zeros(size(xb));
    
    % assemble lists of velocities and coordinates
    X = [X(uselattice); xb(usebdy)];
    Y = [Y(uselattice); zb(usebdy)];
    ux = [ux(uselattice); uxb(usebdy)];
    uy = [uy(uselattice); uyb(usebdy)];


    % interpolate velocities at the measurement points
    Fx = TriScatteredInterp(X,Y,ux);
    Fy = TriScatteredInterp(X,Y,uy);

    uxs(irr) = Fx(xirr,yirr);
    uys(irr) = Fy(xirr,yirr);
end


%% subfunctions

function ind = sub2ind2D(siz,i,j)

% fast conversion of subscripts to matrix indices for 3D arrays

ind = (j-1)*siz(1) + i;