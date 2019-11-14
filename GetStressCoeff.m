function Scoeff = GetStressCoeff(D,Q,c)

% computes on- and off-diagonal coefficients in LBM stress tensor
% D is number of dimensions
% Q is number of velocities
% c is vector of velocity coefficients
% e.g., for D2Q9, D = 2, Q = 9, c = [cx; cy], where cx and cy are the row
% vectors:
% cx  = [ 0, 1, 0, -1, 0, 1, -1, -1, 1]; 
% cy  = [ 0, 0, 1, 0, -1, 1, 1, -1, -1];
% so cx = c(1,:), cy = c(2,:)

Scoeff = zeros(D,D,Q);

for m = 1:D 
    for n = 1:D 
        for i=1:Q             
            Scoeff(m,n,i) = c(m,i) * c(n,i);
        end
    end
end
