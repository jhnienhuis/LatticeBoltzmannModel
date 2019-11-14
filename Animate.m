function Animate(m)


%% load saved data

variables = ['''p'',''iteration'',''time'',''' m.plotquant ''''];
eval([ 'load([m.filedir ''/'' m.filename],' variables ')' ])
eval(['q=' m.plotquant ';']) % q is the quantity to plot

if strcmp(m.plotquant,'velocity') % if we're plotting velocity, calculate its magnitude
    q = sqrt(q.ux.*q.ux + q.uy.*q.uy);
end

m.doTracers = 0;
if isfield(m,'tracerSpacing')
    if m.tracerSpacing % if we're doing tracers
        m.doTracers = 1;
        p.tracerSpacing = m.tracerSpacing;
        load([m.filedir '/' m.filename],'velocity') % we will need x and y components of velocity
    end
end


%% prepare to make movie

% preprocess quantity to be plotted according to boundary conditions
q = PreprocessQuantity(q,p);

% initialize tracer particles and preprocess velocity fields
if m.doTracers
    g = InitializeTracers(p); 

    velocity.ux = PreprocessVelocity(velocity.ux,p);
    velocity.uy = PreprocessVelocity(velocity.uy,p);

end

% determine extreme values over entire run for color plot limits
m.qmin = min(q(:));
m.qmax = max(q(:));

if ~exist([m.filedir '/' m.moviedir],'dir')
    mkdir([m.filedir '/' m.moviedir])
end
 


%% draw and save frames

fig = figure('color','w');

for n = 1:length(iteration)
    
    m.n = n;
    m.it = iteration(n);
            
    if ~mod(m.it/p.saveint,m.frameint) % || n==1 % always save first frame
        SaveFrame(m,p,g,q)
    end
    
    % if we're doing tracers, update tracer locations, even if we aren't
    % saving a frame.
    if m.doTracers
        g = UpdateTracers(p,g,velocity.ux(:,:,n),velocity.uy(:,:,n));
    end

end

close(fig)


%% subfunctions

function SaveFrame(m,p,g,q)

quantity = q(:,:,m.n);
g.bbRegion = find(p.obst);
g.n = m.it;
g.colorlim = [m.qmin m.qmax];

if m.doTracers
    p.doTracers = 1;
end

DrawLBMPlot(p,g,quantity);

print(gcf, ['-d' m.imagetype], [m.filedir '/' m.moviedir '/' m.filename num2str(m.it, '%05d') '.' m.imagetype]);



function q = PreprocessQuantity(q,p)

% force free-slip boundaries to have same value as neighboring row/column
if strcmp(p.bdy.left,'freeslip')
    q(1,:,:) = q(2,:,:);
end

if strcmp(p.bdy.right,'freeslip')
    q(p.lx,:,:) = q(p.lx-1,:,:);
end

if strcmp(p.bdy.lower,'freeslip')
    q(:,1,:) = q(:,2,:);
end

if strcmp(p.bdy.upper,'freeslip')
    q(:,p.ly,:) = q(:,p.ly-1,:);
end


function u = PreprocessVelocity(u,p)

% force obstacles to have zero velocity
u(repmat(p.obst,[1 1 size(u,3)])==1) = 0;

% force free-slip boundaries to have same velocity as neighboring row/column
if strcmp(p.bdy.left,'freeslip')
    u(1,:,:) = u(2,:,:);
end

if strcmp(p.bdy.right,'freeslip')
    u(p.lx,:,:) = u(p.lx-1,:,:);
end

if strcmp(p.bdy.lower,'freeslip')
    u(:,1,:) = u(:,2,:);
end

if strcmp(p.bdy.upper,'freeslip')
    u(:,p.ly,:) = u(:,p.ly-1,:);
end

% pad according to boundary conditions 
u = PadBoundaries(u,p);


function M = PadBoundaries(M,p)

switch p.bdy.left % (actually the first row of the matrix)
    
    case 'noslip'
        M = cat(1,M(1,:,:),M); % pad with extension
        
    case 'freeslip'
        M = cat(1,M(1,:,:),M); % pad with extension
        
    otherwise % boundary will be periodic, the default
        M = cat(1,M(p.lx,:,:),M); 
        
end

switch p.bdy.right % (actually the last row of the matrix)
    
    case 'noslip'
        M = cat(1,M,M(p.lx,:,:)); % pad with extension
        
    case 'freeslip'
        M = cat(1,M,M(p.lx,:,:)); % pad with extension
        
    otherwise % boundary will be periodic, the default
        M = cat(1,M,M(1,:,:)); 
        
end

switch p.bdy.lower % (actually the first column of the matrix)
    
    case 'noslip'
        M = cat(2,M(:,1,:),M); % pad with extension
        
    case 'freeslip'
        M = cat(2,M(:,1,:),M); % pad with extension
        
    otherwise % boundary will be periodic, the default
        M = cat(2,M(:,p.ly,:),M); 
        
end

switch p.bdy.upper % (actually the last column of the matrix)
    
    case 'noslip'
        M = cat(2,M,M(:,p.ly,:)); % pad with extension
        
    case 'freeslip'
        M = cat(2,M,M(:,p.ly,:)); % pad with extension
        
    otherwise % boundary will be periodic, the default
        M = cat(2,M,M(:,1,:)); 
        
end
