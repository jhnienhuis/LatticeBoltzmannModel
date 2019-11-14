function M = FastPadPeriodic(M,n)
% M = FastPadPeriodic(M,n)
% pads a 2D array periodically. n is a 2-element vector containing
% the number of elements to pad with in the 2 directions

% pad rows
M = [M(end-n(1)+1:end,:); M; M(1:n(1),:)];

% pad columns
M = [M(:,end-n(2)+1:end) M M(:,1:n(2))];
