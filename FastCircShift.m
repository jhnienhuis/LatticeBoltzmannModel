function b = FastCircShift(a,p)
%CIRCSHIFT Shift array circularly.
%   B = CIRCSHIFT(A,SHIFTSIZE) circularly shifts the values in the array A
%   by SHIFTSIZE elements. SHIFTSIZE is a vector of integer scalars where
%   the N-th element specifies the shift amount for the N-th dimension of
%   array A. If an element in SHIFTSIZE is positive, the values of A are
%   shifted down (or to the right). If it is negative, the values of A
%   are shifted up (or to the left). 
%
%   Examples:
%      A = [ 1 2 3;4 5 6; 7 8 9];
%      B = circshift(A,1) % circularly shifts first dimension values down by 1.
%      B =     7     8     9
%              1     2     3
%              4     5     6
%      B = circshift(A,[1 -1]) % circularly shifts first dimension values
%                              % down by 1 and second dimension left by 1.
%      B =     8     9     7
%              2     3     1
%              5     6     4
%
%   See also FFTSHIFT, SHIFTDIM, PERMUTE.

%   Copyright 1984-2010 The MathWorks, Inc.  
%   $Revision: 1.11.4.3 $  $Date: 2010/02/25 08:08:46 $
% 
%   Modified by J. T. Perron March 2011 to avoid slow argument checking


sizeA = size(a);
numDimsA = ndims(a);

% Calculate the indices that will convert the input matrix to the desired output
% Initialize the cell array of indices
idx = cell(1, numDimsA);

% Loop through each dimension of the input matrix to calculate shifted indices
for k = 1:numDimsA
    m      = sizeA(k);
    idx{k} = mod((0:m-1)-p(k), m)+1;
end

% Perform the actual conversion by indexing into the input matrix
b = a(idx{:});
