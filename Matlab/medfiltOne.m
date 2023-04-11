function y = medfiltOne(x,n)
%MEDFILT1  One dimensional median filter.
%   Y = MEDFILT1(X,N) returns the output of the order N, one dimensional
%   median filtering of X.  Y is the same size as X; for the edge points,
%   zeros are assumed to the left and right of X.  If X is a matrix,
%   then MEDFILT1 operates along the columns of X.
%
%   If you do not specify N, MEDFILT1 uses a default of N = 3.
%   For N odd, Y(k) is the median of X( k-(N-1)/2 : k+(N-1)/2 ).
%   For N even, Y(k) is the median of X( k-N/2 : k+N/2-1 ).
%

if nargin < 2, n = []; end

% Check if the input arguments are valid
if isempty(n)
  n = 3;
end

nx = length(x);
if rem(n,2)~=1    % n even
    m = n/2;
else
    m = (n-1)/2;  % n odd
end
X = [zeros(m,1); x; zeros(m,1)];
y = zeros(nx,1);

indr = (0:n-1)';
indc = 1:nx;
ind = indc(ones(1,n),1:(nx)) + indr(:,ones(1,nx));
xx = reshape(X(ind),n,nx);
y(1:(nx)) = median(xx,1);
