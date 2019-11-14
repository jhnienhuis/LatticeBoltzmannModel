function omega = GetLESOmega(p,g)
    
% Get a matrix of new values of the relaxation parameter, omega, for the
% Smagorinsky LES turbulence model 
%
% fIn is the array of distributions
% fEq is the array of equilibrium distributions
% D is number of dimensions
% Q is number of velocities
% e.g., for D2Q9, D = 2, Q = 9
% ond and offd are the on- and off-diagonal coefficients in LBM stress
% tensor, computed by GetStressArrays.m
%
% with Chelsea's modifications
%

qadd = 0;
for i = find(g.offd(1,:) ~= 0)
    qadd = qadd + g.offd(1,i)*(g.fNeq(i,:,:));
end
Qo = 2*qadd.*qadd;


for m = 1:p.D % number of diagonal elements
    qadd = 0;
    for i = find(g.ond(m,i) ~= 0) %
        qadd = qadd + g.ond(m,i)*(g.fNeq(i,:,:));
    end
    Qo = Qo + qadd.*qadd;
end


% % ---------------------------- original code --------------------------------------

% % Hou, Sterling, Chen & Doolen (1994) modified Smagorinsky model
% Csqr = C*C;
% S = (sqrt( nu*nu + 18*Csqr*sqrt(Qo) ) - nu) / (6*Csqr); % Theurey thesis, Hou et al
% omega = 1./( 3*( nu + Csqr*S ) + 0.5 );


% % ---------------------------- modified code --------------------------------------

% these lines are needed for both the regular Smagorinsky model & for the Shear-improved model
Qmag= sqrt(Qo);
Csqr= p.C*p.C;
tau0= 3*p.nu*p.A+1/2; % modified for Lax-Wendroff

% if compressible: to strictly account for compressibility effects, should inclue rho
% rho= sum(fIn);

% if incompressible:
rho= 1;

% % Smagorinsky model- e.g., Dong, Sagaut, and Marie 2008
omega= 2./( sqrt(tau0*tau0 + 18*Csqr.*Qmag./rho) + tau0 );  % = 1/tau

% % Shear-improved Smagorinsky model- Jafari & Rahnama 2010, Leveque et al 2007
% % Smag is shear magnitude, calculated somewhere else (I calculate it in Collide.m)
% % I added Smag to the list on inputs

% omega= 2./( sqrt( (tau0 - 3*Csqr*Smag).^2 + 18*Csqr*Qmag./rho ) + tau0 - 3*Csqr*Smag );
