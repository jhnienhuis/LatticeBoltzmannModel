function s = GetStressTensor(omega,nu,fIn,fEq,D,Scoeff)

% omega is the relaxation parameter (= 1/tau)
% nu is kinematic viscosity of LBM
% fIn is the array of distributions
% fEq is the array of equilibrium distributions
% D is number of dimensions (e.g., for D2Q9, D = 2, Q = 9)
% Scoeff is an array of coefficients in LBM stress
% tensor, computed by GetStressCoeff.m
%
% s is the deviatoric stress tensor,
%
% [ sigma_xx - p  sigma_yx     ]
% [ sigma_xy      sigma_yy - p ]
% 
% where sigma_ij are the components of the full stress tensor, 
% and p is the hydrostatic component.
%
% The components of s are:
% s(1,1,:,:) = sxx = sigma_xx - p
% s(2,2,:,:) = syy = sigma_yy - p
% s(2,1,:,:) = sxy = sigma_xy
% s(1,2,:,:) = syx = sigma_yx
% 
% s is symmetric and traceless, such that sxx = -syy and sxy = syx.
%
% The strain rate tensor S = 1/2 * (dua/db + dub/da), with a,b = 1:D, is
% S = s./(2*rho*nu). The components of S are:
%
% S(1,1,:,:) = dux/dx
% S(2,2,:,:) = duy/dy
% S(2,1,:,:) = S(1,2,:,:) = 1/2 * (dux/dy + duy/dx)
%
% (Note that all the above is for D=2 only.)



if D == 2
    [Q lx ly] = size(fIn);
    s = zeros(D,D,lx,ly);
elseif D == 3
    [Q lx ly lz] = size(fIn);
    s = zeros(D,D,lx,ly,lz);
end
    
fdiff = fIn - fEq;
for m = 1:D 
    for n = 1:D
        q = sum( bsxfun( @times, permute(Scoeff(m,n,:),[3 1 2]), fdiff ), 1 );
        s(m,n,:,:) = -3 * nu * omega.*q;
    end
end

