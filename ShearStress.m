function tau = ShearStress(s,ux,uy)

% Rotates the deviatoric stress tensor s at each point into an orientation 
% parallel to the velocity (ux,uy) at that point. Near the boundary, we can 
% use the shear stresses as estimates of the boundary shear stress.
% 
% Returns:
% tau, an lx-by-ly matrix of shear stresses in the same units as s
% nx and ny, lx-by-ly matrices giving the x and y components of unit
% normal vectors to the velocity direction.

% build an array of rotation matrices. To find the shear stress in the
% coordinate parallel to the velocity, we rotate the deviatoric stress
% tensor at each point with a rotation matrix
%
% R = [ cos theta  -sin theta ] = [ ux/u  -uy/u ]
%     [ sin theta   cos theta ]   [ uy/u   ux/u ]
% 
% where u = sqrt(ux^2 + uy^2) is the magnitude of velocity. The new tensor
% is R*S*R', where S is the deviatoric stress tensor at a given point.
%
% For a symmetric (sxy = syx), traceless (sxx = -syy) stress tensor, the
% shear stresses in the new orientation are
%
% sxy' = 1/u^2 * (2*ux*uy*sxx + (ux^2 - uy^2)*sxy)
%
% But note that the assumption of a zero trace only holds for an
% incompressible fluid, whereas the LB method works by allowing some
% compressibility. So there will be small differences between the two
% methods below.
%
% Note also the negative signs on sxx and syy. This is necessary to make
% the computed shear stress in the direction of the velocity match the
% shear stress computed by differencing the velocity in that orientation.


% Method 1: assume symmetry (irrotational) and zero trace (incompressible)

[lx ly] = size(ux);

ux2 = ux.*ux;
uy2 = uy.*uy;

sxx = reshape(s(1,1,:,:),lx,ly);
sxy = reshape(s(2,1,:,:),lx,ly);

tau = (2*ux.*uy.*-sxx + (ux2-uy2).*sxy)./(ux2 + uy2);



% % Method 2: do the full calculation without making assumptions about the
% % stress tensor
% 
% [lx ly] = size(ux);
% 
% u = sqrt(ux.*ux + uy.*uy);
% 
% tau = zeros(size(u));
% 
% for i =1:lx
%    for j = 1:ly
%         R = [ux(i,j) -uy(i,j); uy(i,j) ux(i,j)]/u(i,j);
%         S = [-s(1,1,i,j) s(1,2,i,j); s(2,1,i,j) -s(2,2,i,j)]; % [-sxx syx; sxy -syy]
%         Sprime = R*S*R';
%         tau(i,j) = Sprime(2,1);
%    end
% end
