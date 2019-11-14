function LBM(varargin)
%
%changed: mar17, 2013 by Jaap
% added Bouzidi to possible boundary conditions 
%changed: mar31, 2013 by Jaap
% subtle stability/reality checks
%
%future:
% model in subversion
% netcdf
% compare bouzidi shear stresses with original shear stress
% clean up MicroBC

%get timing of run
maxNumCompThreads(8 ); tic;

%parse input
if nargin && isstruct(varargin{1}),
    p = varargin{1};
elseif nargin && ischar(varargin{1}),
    p = struct('savename',varargin{1});
else
    p = struct;
end
%create p struct (all parameters that not change during run)
[p]   = InitializeStruct(p);  

fprintf('Initializing: %s \n', p.savename)

%create g struct (all parameters that change during run)
[p,g] = InitializeLBM(p);

fprintf('Running: %s \n', p.savename)

%run LB model
[p,g] = RunLBM(p,g);

fprintf('Finalizing: %s \n', p.savename)


p.RunTime = toc;

%write and plot commands
FinalizeLBM(p,g);


end