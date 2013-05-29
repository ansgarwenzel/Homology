function H = householder(b, k)
% H = householder(b, k)
% Atkinson, Section 9.3, p. 611
% b is a column vector, k an index < length(b)
% Constructs a matrix H that annihilates entries
% in the product H*b below index k

% $Id: householder.m,v 1.1 2008-01-16 15:33:30 mike Exp $
% M. M. Sussman

n = length(b);
d(:,1) = b(k:n);
if d(1)>=0
  alpha = -norm(d);
else
  alpha =  norm(d);
end

if alpha==0
  H=eye(n);
  return
end

lenD=length(d);
v=zeros(lenD,1);

v(1,1)=sqrt(.5*(1-d(1)/alpha));
p=-alpha*v(1);
v(2:lenD)=d(2:lenD)/(2*p);
w=[zeros(k-1,1);v];
H=eye(n)-2*w*w';