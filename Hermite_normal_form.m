function [H U] = Hermite_normal_form(A)
% (column-reduced) Hermite normal form over integers
%
% Input: Matrix A with integers as entries
% Output: H is the lower-triangular Hermite normal form of A obtained by
% the following three types of column operations
% (i) exchange two columns
% (ii) multiply a column by -1
% (iii) Add an integral multiple of a column to another column
%
% U is the unitary matrix such that AU = H
%
%
% A full-rank square matrix with integer entries is said to be
% in Hermite normal form if 
% (i) it is lower triangular, 
% (ii) the diagonal entries are positive and 
% (iii) all entries to the left of a diagonal entries are non-negative
% and strictly less than the corresponding diagonal entry.
%
% In general, an m x n integral matrix H of rank r is in Hermite normal form
% if
% (i) the first r columns are nonzero. The n-r columns on the right are zero.
% (ii) for i=1,2,...,r, we have r pivots p(1), p(2), ..., p(r), such that
% p(1) < p(2) < ...< P(r) and H(p(i),i) is the first nonzero
% entry in column i.
% (iii) for i=1,2,...,r, H(p(i), i) > 0.
% (iv) for i=1,2,...,r, and j=1,...,i, H(p(i), i) > H(p(i),j) >=0.
%
% For example, the matrix 
%
% A =
% 
% 9 6 0 -8 0
% -5 -8 0 0 0
% 0 0 0 4 0
% 0 0 0 -5 0
% 
% can be column-reduced to the following HNF
%
% H =
% 
% 1 0 0 0 0
% 1 2 0 0 0
% 28 36 84 0 0
% -35 -45 -105 0 0
%
%
[m n] = size(A);
% Input matrix A has m rows and n columns
A = round(A); % make sure that the input is integer-valued
H = A; % initialize H to be A
U = eye(n,n); % initialize U to the identity matrix
r = min(m,n);
p = 0; % the location of pivot
for i = 1:r;
  S = find((H(i,(p+1):n)) ~= 0); % Find the nonzero entries in row i
  if ~isempty(S)
    p = p+1;
    FINISHED = 0;
    while ~FINISHED
      S = p-1 + find((H(i,p:n)) ~= 0);
      [dummy k] = min(abs(H(i,S))); % find the smallest non-zero entry
      if S(k) ~= p
        index = 1:n; index(p) = S(k); index(S(k)) = p;
        H = H(:,index); % exchange columns i and k
        U = U(:,index); % exchange columns i and k
      end
        for j = (p+1):n
        q = round(H(i,j)/H(i,p));
        H(:,j) = H(:,j) - q*H(:,p);
        U(:,j) = U(:,j) - q*U(:,p);
      end
      if isempty(find( H(i,(p+1):end) ~= 0) )
        FINISHED = 1;
        if H(i,p) < 0 % flip the sign of H(i,;) if necessary
          H(:,p) = -H(:,p);
          U(:,p) = -U(:,p);
        end
 % reduce the entries to the left of H(i,i)
        for j = 1:(p-1)
          q = floor(H(i,j)/H(i,p));
          H(:,j) = H(:,j) - q*H(:,p);
          U(:,j) = U(:,j) - q*U(:,p);
        end
 end % if
 end % while
 end % if
end % for