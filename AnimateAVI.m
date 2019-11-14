function AnimateAVI(out)

do_curl = 0;
do_save = 1;
do_velocity = 1;

if do_save,
    % load saved data
    if isa(out,'char'),
        %assume out is a .mat file located on the file path
        out = load(out);
    elseif isa(out,'struct'),
        if any(~ismember({'iteration','time','p','velocity'},fieldnames(out))),
            error('Loaded struct does not contain all neccessary fields.')
        end
    end
    %out.p.savedir = [dropbox filesep 'WaveRipple' filesep 'data_lb'];
    out.p.savedir = cd;
    
    if ~exist(out.p.savedir,'dir'),
        mkdir(out.p.savedir)
    end
    
    %save([out.p.savedir filesep out.p.savename '.mat'],'-struct','out','p','bedshearstress','iteration','time')
end

if do_velocity,
    q = sqrt(out.velocity.ux.*out.velocity.ux + out.velocity.uy.*out.velocity.uy);
    
    % determine extreme values over entire run for color plot limits

out.p.qmax = max(q(:));
out.p.qmin = 0;
elseif do_curl,
    q = zeros(size(out.velocity.ux));
    for i=1:size(out.velocity.ux,3),
        q(:,:,i) = curl(out.p.x0*(1:out.p.ly),out.p.x0*(1:out.p.lx),out.velocity.ux(:,:,i),out.velocity.uy(:,:,i));
    end
    
    % determine extreme values over entire run for color plot limits
out.p.qmax = 2*std(q(:));
out.p.qmin = -out.p.qmax;
end
%% prepare to make movie

% preprocess quantity to be plotted according to boundary conditions
q = PreprocessQuantity(q,out.p);

% initialize tracer particles and preprocess velocity fields
out.p.tracerSpacing = 10;
g = InitializeTracers(out.p); 

out.velocity.ux = PreprocessVelocity(out.velocity.ux,out.p);
out.velocity.uy = PreprocessVelocity(out.velocity.uy,out.p);


 


%% draw and save frames

if do_curl,
    writerOBJ = VideoWriter([out.p.savedir filesep out.p.savename 'curl.avi']);
elseif do_velocity,
    out.p.savedir = cd;
    %out.p.savedir = 'D:\Users\jaap\Dropbox (MIT)\';
    %out.p.savedir = [dropbox filesep 'FieldRipples' filesep 'data_lb'];
    %out.p.savedir = 'D:\FieldRipples\data_lb';
    writerOBJ = VideoWriter([out.p.savedir filesep out.p.savename '.avi']);
end
writerOBJ.FrameRate = 15;
open(writerOBJ)

fig = figure('color','w','visible','off');

for n = 1:numel(out.time);
    
    out.p.n = n;
    out.p.it = out.iteration(n);
            
    SaveFrame(out,g,q);
    
    writeVideo(writerOBJ,imcapture(fig))
        
    disp(['Writing Frame # ' num2str(out.p.it)])

    % if we're doing tracers, update tracer locations, even if we aren't
    % saving a frame.
    g = UpdateTracers(out.p,g,out.velocity.ux(:,:,n),out.velocity.uy(:,:,n));
    


end

close(fig)
close(writerOBJ)


%% subfunctions

function SaveFrame(out,g,q)

quantity = q(:,:,out.p.n);
g.bbRegion = find(out.p.obst);
g.n = out.p.it;
g.colorlim = [out.p.qmin out.p.qmax];

out.p.doTracers = 1;

DrawLBMPlot(out.p,g,quantity);

end

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
end

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
end



end