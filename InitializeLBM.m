function [p,g] = InitializeLBM(p)
%fill G structure, all LB specific/derived parameters

%% scaling conversions and lattice constants

p.nu    = 1/(p.rho0/p.mu0*p.x0^2/p.t0); % LB kinematic viscosity from scaling parameters
%p.nu    = p.t0/(p.x0^2)*(1/p.Re);

g.omega = 1. / (3*p.nu*p.A+1./2.); % relaxation parameter, starting value (=1/tau)

% scalings to LBM units
p.m0 = p.rho0*p.x0^3;          % mass scale
p.p0 = p.m0/p.x0/p.t0^2;       % stress scale 
p.f0 = p.rho0*p.x0/p.t0^2;     % force scale
p.u0 = p.x0/p.t0;              % velocity scale


g.n = 0; % iteration count

% g.taub = zeros(p.lx,1);

%   7  3  6  
%    \ | /   
%   4- 1 -2  
%    / | \   
%   8  5  9  

% D2Q9 LATTICE CONSTANTS 
g.t   = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36]; 
g.cx  = [ 0, 1, 0, -1, 0, 1, -1, -1, 1]; 
g.cy  = [ 0, 0, 1, 0, -1, 1, 1, -1, -1]; 
g.opp = [ 1, 4, 5, 2,  3, 8, 9,  6,  7];


%% external forces

if length(p.F) == 1 % if F is a scalar, assumed to be constant force in the x direction
    p.F = p.F*[ones(p.N+1,1) zeros(p.N+1,1)];
elseif length(p.F) == 2 % constant force in [x y] directions
    p.F = repmat(p.F(:)',[p.N+1,1]);
elseif length(p.F) == p.N+1 && size(p.F,2)==1 % time series of forces in x direction
    p.F = [p.F(:) zeros(p.N+1,1)];
else % time series of forces, with each row being [x y] directions
    p.F = p.F;
end

%% stress arrays and turbulence model

if isfield(p,'doLES') && p.doLES
    [g.ond,g.offd] = GetStressArrays(p.D,p.Q,[g.cx; g.cy]);
end

g.Scoeff = GetStressCoeff(p.D,p.Q,[g.cx; g.cy]);


if isfield(p,'MRT') && p.MRT
    
    % transformation matrix for two relaxation time
    % r = (rho, e, epsilon, jx, qx, jy, qy, pxx, pxy)
    % see Lallemand and Luo, JCP 2003
    
    g.M = [  1  1  1  1  1  1  1  1  1;
            -4 -1 -1 -1 -1  2  2  2  2;
             4 -2 -2 -2 -2  1  1  1  1;
             0  1  0 -1  0  1 -1 -1  1;
             0 -2  0  2  0  1 -1 -1  1;
             0  0  1  0 -1  1  1 -1 -1;
             0  0 -2  0  2  1  1 -1 -1;
             0  1 -1  1 -1  0  0  0  0;
             0  0  0  0  0  1 -1  1 -1];
    
    if isfield(p,'s2')
        g.s2 = p.s2; % bulk viscosity
        g.s3 = p.s3;
        g.s57 = p.s57;  % note! 5 and 7 ought to be same value?
    else
        % Lallemand and Luo JCP 2003 values
        g.s2 = 1.5; % bulk viscosity
        g.s3 = 1.4;
        g.s57 = 1.5;
    end
    
    % Lallemand and Luo PRE 2000, fig 1 values
    % g.s2 = 1.64; % FIXME bulk viscosity
    % g.s3 = 1.54;
    % g.s5 = 1.9;
    % g.s7 = 1.9;
    
end


%% boundary conditions

% obstacles, no-slip or bouzidi conditions are applied
p.obst  = zeros(p.lx,p.ly); 
p.obst(p.solid == 1) = 1;

% for free-slip BCs

%   7  3  6  
%    \ | /   
%   4- 1 -2  
%    / | \   
%   8  5  9  

% default is periodic
g.fsN = 0;
g.fsS = 0;
g.fsE = 0;
g.fsW = 0;

switch p.bdy.left
    
    case 'noslip'
        p.obst(1,:) = 1;
        
    case 'bouzidi'
        p.obst(1,:) = 1;
        
    case 'freeslip'
        g.fsW  = 1;
        
        g.fsWind = [1, 4, 3, 2, 5, 7, 6, 9, 8];
        
        slip = zeros(p.lx,p.ly);
        slip(1,1:p.ly) = 1;
        g.fsWrecip = find(slip);

        slip = zeros(p.lx,p.ly);
        slip(2,1:p.ly) = 1;
        g.fsWdonor = find(slip);
        
    otherwise
        % boundary will be periodic, the default
        
end

switch p.bdy.right
    
    case 'noslip'
        p.obst(p.lx,:) = 1;
        
    case 'bouzidi'
        p.obst(p.lx,:) = 1;
        
    case 'freeslip'
        g.fsE  = 1;
        
        g.fsEind = [1, 4, 3, 2, 5, 7, 6, 9, 8];

        slip = zeros(p.lx,p.ly);
        slip(p.lx,1:p.ly) = 1;
        g.fsErecip = find(slip);

        slip = zeros(p.lx,p.ly);
        slip(p.lx-1,1:p.ly) = 1;
        g.fsEdonor = find(slip);
        
    otherwise
        % boundary will be periodic, the default
        
end

switch p.bdy.lower
    
    case 'noslip'
        p.obst(:,1) = 1;
        
    case 'bouzidi'
        p.obst(:,1) = 1;
        
    case 'freeslip'
        g.fsS  = 1;
        
        g.fsSind = [1, 2, 5, 4, 3, 9, 8, 7, 6];
        
        slip = zeros(p.lx,p.ly);
        slip(1:p.lx,1) = 1;
        g.fsSrecip = find(slip);

        slip = zeros(p.lx,p.ly);
        slip(1:p.lx,2) = 1;
        g.fsSdonor = find(slip);
        
    otherwise
        % boundary will be periodic, the default
        
end

switch p.bdy.upper
    
    case 'noslip'
        p.obst(:,p.ly) = 1;
        
    case 'bouzidi'
        p.obst(:,p.ly) = 1;
        
    case 'freeslip'
        g.fsN  = 1;
        
        g.fsNind = [1, 2, 5, 4, 3, 9, 8, 7, 6];
        
        slip = zeros(p.lx,p.ly);
        slip(1:p.lx,p.ly) = 1;
        g.fsNrecip = find(slip);

        slip = zeros(p.lx,p.ly);
        slip(1:p.lx,p.ly-1) = 1;
        g.fsNdonor = find(slip);
        
    otherwise
        % boundary will be periodic, the default
        
end


% create list of matrix indices for which no-slip boundary condition applies
g.bbRegion = find(p.obst); 


%% initial conditions and preallocation of arrays

% Zero velocity (rho=0, u=0) ==> fIn(i) = t(i) 
g.fzero = reshape( g.t' * ones(1,p.lx*p.ly), p.Q, p.lx, p.ly); 
g.fIn = reshape( g.t' * ones(1,p.lx*p.ly), p.Q, p.lx, p.ly); 

% % constant velocity
% ux = p.U * ones(p.lx,p.ly);
% uy = zeros(p.lx,p.ly);
% g.rho = 1;
% for i=1:p.Q
%     cu = 3*(g.cx(i)*ux+g.cy(i)*uy);
%     g.fIn(i,:,:) = g.rho .* g.t(i) .* ...
%                    ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2) );
% end
% ux(g.fsNRegion) = 0;   % Chelsea
% uy(g.fsNRegion) = 0;

% calculate initial density and velocity
g = GetRhoU(p,g);

% Preallocate arrays
g.fEq = zeros(p.Q,p.lx,p.ly);
g.fOut = zeros(p.Q,p.lx,p.ly);
g.fNeq = zeros(p.Q,p.lx,p.ly);
%g.force = zeros(p.Q,p.lx,p.ly);
g.force = zeros(p.Q,1);  % uniform force

if isfield(p,'MRT') && p.MRT
    g.rIn = zeros( p.Q, p.lx*p.ly );
    g.rEq = zeros( p.Q, p.lx*p.ly );
    g.rNeq = zeros( p.Q, p.lx*p.ly );
    g.rOut = zeros( p.Q, p.lx*p.ly );
end

% calculate initial stresses
g = GetStress(p,g);

% /J Initialize boundary conditions based on initial bed elevations
% Changed since 3/17/13. 
if any(strcmpi(struct2cell(p.bdy),'bouzidi')),
    g = GetBouzidi(p,g);
end

% get initial body force
g = GetForce(p,g);


% create arrays for saved output, if applicable, and save initial conditions

if isfield(p,'saveint') && p.saveint % If we're saving results
    
    lsave = p.N/p.saveint + 1; % number of frames we will save. The +1 is because we also save the initial state
    
    g.iteration = zeros(lsave,1);
    g.time = zeros(lsave,1);
    
    if any(ismember(p.savequant,'velocity'))
        
        g.velocity.ux = zeros(p.lx,p.ly,lsave);
        g.velocity.uy = zeros(p.lx,p.ly,lsave);
        g.velocity.ux(:,:,1) = g.ux*p.u0; % x velocity in m/s
        g.velocity.uy(:,:,1) = g.uy*p.u0; % y velocity in m/s
        
    end
    
    if any(ismember(p.savequant,'shearstress'))
        
        g.shearstress = zeros(p.lx,p.ly,lsave);
        g.shearstress(:,:,1) = g.tau*p.p0; % shear stress in Pa
        
    end
    
    if any(ismember(p.savequant,'pressure'))
        
        g.pressure = zeros(p.lx,p.ly,lsave);
        g.pressure(:,:,1) = (1/3)*g.rho*p.p0; % pressure in Pa
        
    end
    
    if any(ismember(p.savequant,'bedshearstress'))
        
        g.bedshearstress = zeros(p.lx,1,lsave);
        
    end
    
end
