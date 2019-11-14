function g = GetBouzidi(p,g)

% applies modified bounce-back BC described by Bouzidi et al. (2001) for a
% boundary with known location

%% find all lattice links that cross no-slip boundary

solid = zeros(p.lx,p.ly);
solid(g.bbRegion) = 1; 
fluid = ~solid; % ignores free-slip boundaries
% fluid = padarray(fluid,[1 0],'circular'); % periodic x boundaries
fluid = FastPadPeriodic(fluid,[1 0]); % periodic x boundaries
fluid = bwperim(fluid,8); % fluid nodes that have at least one 8-connected neighbor across the boundary
fluid = fluid(2:end-1,:); % recover original x dimension
fluid(:,p.ly) = 0; % exclude upper y boundary

% generate a list of links with x and y coordinates of fluid node and the
% direction of the link to the solid node

indf = []; % matrix indices of fluid nodes
% inds = []; % matrix indices of solid nodes
dir = []; % directions of links from fluid to solid nodes

for i=1:p.Q 
    newindf = find(fluid & FastCircShift(solid, [-g.cx(i),-g.cy(i)])); 
%     newinds = find(solid & FastCircShift(fluid, [ g.cx(i), g.cy(i)]));
    indf = cat(1,indf,newindf);
%     inds = cat(1,inds,newinds);
    dir = cat(1,dir,i*ones(size(newindf)));
end

[y,x] = meshgrid(1:p.ly,1:p.lx);
xf = x(indf); % x coordinates of fluid nodes
yf = y(indf); % y coordinates of fluid nodes

% xs = x(inds); % x coordinates of solid nodes
% ys = y(inds); % x coordinates of solid nodes
xs = xf + g.cx(dir)';
ys = yf + g.cy(dir)';
xs(xs<1) = xs(xs<1) + p.lx; % apply periodic x BC
xs(xs>p.lx) = xs(xs>p.lx) - p.lx; % apply periodic x BC
% inds = sub2ind([p.lx p.ly],xs,ys); % matrix indices

zbf = p.bed(xf); % elevation of boundary at x=xf
zbs = p.bed(xs); % elevation of boundary at x=xs

% upstream fluid node (neighbor of main fluid node in direction away from the boundary)
xf2 = xf - g.cx(dir)';
yf2 = yf - g.cy(dir)';
xf2(xf2<1) = xf2(xf2<1) + p.lx; % apply periodic x BC
xf2(xf2>p.lx) = xf2(xf2>p.lx) - p.lx; % apply periodic x BC
% f2ind = sub2ind([p.lx p.ly],xf2,yf2); % matrix indices


% calculate linear indices to be used in ApplyBC.m
siz = [p.Q p.lx p.ly];

ns.fd = sub2ind3D(siz,dir,xf,yf); % distributions at fluid node in the direction of the boundary
ns.fu = sub2ind3D(siz,g.opp(dir)',xf,yf); % distributions at fluid node in the direction away from the boundary
ns.f2d = sub2ind3D(siz,dir,xf2,yf2); % distributions at upstream fluid node in the direction of the boundary
ns.fout = sub2ind3D(siz,g.opp(dir)',xs,ys); % outgoing directions of solid nodes

%% find Bouzidi's q for each link

ns.q = abs(yf - zbf) ./ (abs(yf - zbf) + abs(ys - zbs));

% for vertical links (fluid and solid node have same x coordinate), it's
% much simpler:
vertical = xf == xs;
ns.q(vertical) = yf(vertical) - zbf(vertical);

% % check q values by plotting intersections of links with boundary
% xb = xf + g.cx(ns.fout)'.*ns.q;
% yb = yf + g.cy(ns.fout)'.*ns.q;
% 
% figure
% imagesc(fluid); axis image; hold on
% plot(g.z,1:p.lx,'-w')
% for i=1:length(xf)
% %     plot([yf(i) ys(i)],[xf(i) xs(i)],'-ow')
%     plot([yf2(i) ys(i)],[xf2(i) xs(i)],'-w')
% end
% plot(yb,xb,'sy')

%% assign output

g.ns = ns;


%% subfunctions

function ind = sub2ind3D(siz,i,j,k)

% fast conversion of subscripts to matrix indices for 3D arrays

ind = (k-1)*siz(1)*siz(2) + (j-1)*siz(1) + i;