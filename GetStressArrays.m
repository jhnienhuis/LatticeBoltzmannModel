function [ond,offd] = GetStressArrays(D,Q,c)

% computes on- and off-diagonal coefficients in LBM stress tensor
% D is number of dimensions
% Q is number of velocities
% c is vector of velocity coefficients
% e.g., for D2Q9, D = 2, Q = 9, c = [cx; cy], where cx and cy are the row
% vectors:
% cx  = [ 0, 1, 0, -1, 0, 1, -1, -1, 1]; 
% cy  = [ 0, 0, 1, 0, -1, 1, 1, -1, -1];
% so cx = c(1,:), cy = c(2,:)

ond = zeros(D,Q);
offd = zeros(D,Q);

odm = 1; % off-diagonal m index

for m = 1:D 
    for n = 1:D 
        for i=1:Q 
            
            em = c(m,i);
            en = c(n,i);
            coeff = em*en;
            if m == n
                ond(m,i) = coeff;
            elseif m > n
                offd(odm,i) = coeff;
            end
        end

        if m == n
            % do nothing
        elseif m > n 
            odm = odm+1; 
        end
    end
end

% for D2Q9
% ond =  [0     1     0     1     0     1     1     1     1; ...
%         0     0     1     0     1     1     1     1     1];
% offd=  [0     0     0     0     0     1    -1     1    -1; ...
%         0     0     0     0     0     0     0     0     0];
