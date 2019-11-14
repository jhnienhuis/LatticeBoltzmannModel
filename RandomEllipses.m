function obst = RandomEllipses(p,xspacing,yspacing,randomness,r,fracsdr,aspect,mask,periodic)

if nargin < 7
    aspect = 1; % if aspect >1, ellipses will have long axes in x direction
end

if nargin < 8
    mask = zeros(p.lx,p.ly); % no ellipse centers where mask == 1
end

if nargin < 9
    periodic = 1; % will ellipses be periodic over the domain?
end


% evenly spaced locations
[yc,xc] = meshgrid(1:yspacing:p.ly,1:xspacing:p.lx);

xc = xc(:);
yc = yc(:);

% add random component to locations
yc = yc + randomness*yspacing*(rand(size(yc))-0.5);
xc = xc + randomness*xspacing*(rand(size(xc))-0.5);

% periodic boundaries
xc = mod( xc - 0.5, p.lx ) + 0.5;
yc = mod( yc - 0.5, p.ly ) + 0.5;


% Eliminate center points whose closest lattice site has mask==1

% round coordinates to the nearest lattice site
xcr = round(xc);
ycr = round(yc);

% convert to linear indices
icr = sub2ind([p.lx p.ly],xcr,ycr);

xc = xc(~mask(icr));
yc = yc(~mask(icr));


% mean radii
r = r*ones(size(xc));
dr = r.*fracsdr.*randn(size(r)); % random deviations from mean r with a fractional standard deviation of fracsdr
r = r + dr;

% semimajor and minor axes
b = 2*r/(aspect+1);
a = b*aspect;


if periodic
    % tile x and y matrix coordinates 3x3 for periodic boundaries
    x = 1:p.lx;
    y = 1:p.ly;
    [Y,X] = meshgrid([y-p.ly y y+p.ly],[x-p.lx x x+p.lx]);

    obst = zeros(size(X));

    for i = 1:length(r)
        d = ((X-xc(i))/a(i)).^2 + ((Y-yc(i))/b(i)).^2; % ellipse function
        obst(d<=1) = 1;
    end

    obst = mat2cell(obst, p.lx*[1 1 1], p.ly*[1 1 1]);
    temp = zeros(p.lx,p.ly);
    for i=1:9
        temp = temp + obst{i};
    end
    
    obst = temp > 0;

else
    
    [Y,X] = meshgrid(1:p.ly,1:p.lx);

    obst = zeros(size(X));

    for i = 1:length(r)
        d = ((X-xc(i))/a(i)).^2 + ((Y-yc(i))/b(i)).^2; % ellipse function
        obst(d<=1) = 1;
    end
    
end

