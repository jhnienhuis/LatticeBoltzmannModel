function [p] = InitializeStruct(q)
%fill P structure, all fluid/flow/simulation parameters


%% Parameters in Lattice Boltzmann Model

p.D = 2;                        % number of dimensions (p.D=2 for 2-D flow)
p.Q = 9;                        % number of directions (typically p.Q=9 for 2-D flow)

p.doLES = 1;                    % 1 to use Smagorinsky large-eddy turbulence model, 0 not
p.C = 0.18;                     % corrected Smagorinsky coefficient

p.A = 0.999;                    % Lax-Wendroff parameter for streaming step. Set to 1 to do standard finite difference

p.MRT = 1;                      % 1 to do multiple relaxation times for high Re flow, 0 not
p.s2 = 1.1;                     % used in MRT model
p.s3 = 1.1;                     % used in MRT model
p.s57 = 1.1;                    % used in MRT model


p.rho0 = 998.2;                 % fluid density (kg/m^3) (water at STP is 998.2)
p.mu0 = 0.00100;                % fluid dynamic viscosity (Pa s) (water at STP is 0.001)

p.x0 = 5e-03;                   % grid spacing (meters per lattice unit)
p.t0 = 1e-04;                   % time step (seconds per LBM step)

p.T = 1;                       % desired duration of run (seconds)
p.N = ceil(p.T/p.t0);           % number of LBM iterations to do

%% Filenames and output

p.saveint = 1000;                % save solution to memory every saveint iterations. To save at a set time interval timeint, use p.saveint = round(timeint/p.t0). To not save, use p.saveint = 0
p.writeint = 1000;             % write cumulative results to disk every writeint iterations. To write to disk only at end, use p.writeint = 0
p.savedir = ['D:' filesep 'data_lb']; % directory to save results in
p.savename = 'FieldRipple25';          % filename for .mat file
p.savequant = {'velocity','bedshearstress'}; % quantities to save. Options are 'velocity','shearstress','pressure','bedshearstress' (last one only saves 1-D data)

p.plotint = 0;               % plot solution every plotint iterations. To not plot, use p.plotint = 0
p.plotquant = 'velocity';      % quantity to plot. Options are 'velocity','shearstress','pressure'
p.plotVectors = 1;             % plot velocity vectors
p.Vectorspacing = 4;           % velocity vector spacing in lattice units
p.plotStreamlines = 0;         % plot flow streamlines

%% Domain and boundary conditions
p.lx = 400;                    % number of points in x direction
p.ly = 100;                    % number of points in y direction

p.bdy.left  = 'periodic';      % left boundary condition. Options are 'periodic', 'noslip', 'freeslip' or 'bouzidi' (interpolated no slip).
p.bdy.right = 'periodic';      % right boundary condition.
p.bdy.upper = 'freeslip';      % upper boundary condition.
p.bdy.lower = 'bouzidi';        % lower boundary condition.

p.bed = []; % which bed elevation model to load?
%% Error checking and input finalization

q_names = fieldnames(q);

for i = 1:numel(q_names)
    %fill q field in p field
    p.(q_names{i}) = q.(q_names{i});
    
    %change savename such that no overwriting
    p.savename = [p.savename '_' q_names{i} num2str(q.(q_names{i}))];
    
end

%% Dependent Variables

p.u0 = 0.6;
p.beta = -7.5* (p.u0/p.lx);                  %factor for horizontal pressure gradient, negative gradient is positive velocity (see Zhang & Kwok Phys Rev E 2006)

%no body force for unidirectional flow
p.F = 0;

%number of calculations
p.NCalc = p.ly*p.lx*p.N;
%% Solid Bed Geometry

BED = 'currentripple';

switch BED,
    case 'none'
        p.bed = zeros(p.lx,1);
        
    case 'flat'
        zmin = 1;
        p.bed = zmin*ones(p.lx,1);
        
        
    case 'sinusoid'
        %a trochoidal ripple shape with a wavelength equal to the length of the domain divided by ncyc
        
        ncyc = 1;                       % number of cycles in the x direction
        a = 0.15;                       % bedform aspect ratio (height as a fraction of wavelength)
        zmin = 1;                       % minimum elevation
        
        p.bed = zmin + (p.lx-1)*a*(1 - abs(sin(2*pi*((1:p.lx)-p.lx/2)'/(p.lx/(ncyc/2))))); % sharp crests, smooth troughs
        
        
    case 'currentripple'
        
        
        ncyc = 3;
        l = 100; 
        
        eta = 0.1*l;
        lee_slope = 30;
        
        bed = zeros(1,l);
        
        
        l_lee = round(eta./tan(deg2rad(lee_slope)));
        
        bed(1:l_lee) = eta*linspace(1,0,l_lee);
        
        bed((l_lee+1):l) = eta*0.5*(1+tanh(4*linspace(-0.3,0.3,l-l_lee))/tanh(1.2));
        
        p.bed = (1+circshift(bed,[1 -l_lee]))';
        
        %make sure p.bed(1) is 1
        p.bed = p.bed + (ceil(p.bed(1))-p.bed(1));
        
        p.bed = [p.bed; ones(p.wlength,1); p.bed; ones(p.lx-200-p.wlength,1)];
        
end


[Y,~] = meshgrid(1:p.ly,1:p.lx);
p.solid = zeros(p.lx,p.ly);     % 1 where there is a solid object in the domain, zero where there is fluid
p.solid(Y<repmat(p.bed(:),[1 p.ly])) = 1;
end