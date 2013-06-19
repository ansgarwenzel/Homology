function [S] = smith(A)

% Smith normal form of an integer matrix.
%
% [U,S,V] = smith(A) returns integer matrices U, S, and V such that
% A = U*S*V',
% S is diagonal and nonnegative, S(i,i) divides S(i+1,i+1) for all i,
% det U =+-1, and det V =+-1.
% s = smith(A) just returns diag(S).
% Uses function ehermite.
% [U,S,V] = smith(A);
%
% This function is in some ways analogous to SVD.

% John Gilbert, 415-812-4487, December 1993
% gilbert@parc.xerox.com 
% Xerox Palo Alto Research Center 

if round(A) ~= A
error('Requires integer input.');
end

% This looks much like an SVD algorithm that first bidiagonalizes
% A by Givens rotations and then chases zeros, except for
% the construction of the 2 by 2 elementary transformation.

[m,n] = size(A);
S = A;
%U = eye(m);
%V = eye(n);

% Bidiagonalize S with elementary Hermite transforms.

% [U,S,V] = bidiag(S);
% U = round(U);
% S = round(S);
% V = round(V);

count_zero_columns = 0;
count_zero_rows = 0;

for i = 0:(n-1),
    if (sum(S(:,(n-i))==0)==m),
        S(:,(n-i):n) = circshift(S(:,(n-i):n), [0,-1]);
        count_zero_columns = count_zero_columns + 1;
    end
end

for i = 0:(m-1),
    if (sum(S((m-i),:)==0)==n),
        S((m-i):m,:) = circshift(S((m-i):m,:), [-1,0]);
        count_zero_rows = count_zero_rows + 1;
    end
end

H = S(1:(m-count_zero_rows),1:(n-count_zero_columns));
[o,p] = size(H);

for j = 1:p,
    if H(j,j)==0,
        y = find(H((j+1):o,j),1);
        H = [H(1:(j-1),:);H((j+y),:);H(j:(y-1+j),:);H((y+1+j):o,:)];
    end   
    for i = (j+1):o,
        if H(i,j),
            x = gcd(H(i,j),H(j,j));
            H(i,j:p) = (H(j,j)/x)*H(i,j:p) - (H(i,j)/x)*H(j,j:p);
        end
    end
end

%make diagonal direcly. needs to be changed somewhat.
%for testdummy = 2:p,
for i = 1:p-1,
    for j = i+1:p,
        if H(i,j),
            x = gcd(H(i,j),H(i,i));
            H(:,j) = (H(i,i)/x)*H(:,j) - (H(i,j)/x)*H(:,i);
        end
    end
end
%end

S(1:o,1:p) = H;

% for j = 1:min(m,n),
% % Zero column j below the diagonal.
% for i = j+1:m
%    % i = m+j+1-k;
%    S = round(S);
% if S(i,j)
% % Construct an elementary Hermite transformation E
% % to zero S(i,j) by combining rows i and j.
% E = ehermite(S(j,j),S(i,j));
% % Apply the transform to S and U.
% S([j i],:) = E * S([j i],:);
% %U(:,[j i]) = U(:,[j i]) / E;
% end;
% end;
% % Zero row j after the superdiagonal.
% for i = j+2:n
% if S(j,i)
% % Construct an elementary Hermite transformation E
% % to zero S(j,i) by combining columns j+1 and i.
% E = ehermite(S(j,j+1),S(j,i));
% % Apply the transform to S and V.
% S(:,[j+1 i]) = S(:,[j+1 i]) * E';
% %V(:,[j+1 i]) = V(:,[j+1 i]) / E;
% end;
% end;
% end;

% Now S is upper bidiagonal.
% Chase the superdiagonal nonzeros away.
% D = diag(S,1);
% while any(D)
% b = min(find(D));
% % Start chasing bulge at first nonzero superdiagonal element.
% 
% % To guarantee reduction in S(b,b), first make S(b,b) positive
% % and make S(b,b+1) nonnegative and less than S(b,b).
% if S(b,b) < 0
% S(b,:) = -S(b,:);
% %U(:,b) = -U(:,b);
% end;
% q = floor(S(b,b+1)/S(b,b));
% E = [1 0 ; -q 1];
% S(:,[b b+1]) = S(:,[b b+1]) * E';
% %V(:,[b b+1]) = V(:,[b b+1]) / E;
% 
% if S(b,b+1)

% Zero the first nonzero superdiagonal element
% using columns b and b+1, to start the bulge at S(b+1,b).
% E = ehermite(S(b,b),S(b,b+1));
% S(:,[b b+1]) = S(:,[b b+1]) * E';
% %V(:,[b b+1]) = V(:,[b b+1]) / E;
% for j = 1:min(m,n)
% if j+1 <= m
% % Zero S(j+1,j) using rows j and j+1.
% E = ehermite(S(j,j),S(j+1,j));
% S([j j+1],:) = E * S([j j+1],:);
% %U(:,[j j+1]) = U(:,[j j+1]) / E;
% end
% if j+2 <= n
% % Zero S(j,j+2) using columns j+1 and j+2.
% E = ehermite(S(j,j+1),S(j,j+2));
% S(:,[j+1 j+2]) = S(:,[j+1 j+2]) * E';
% %V(:,[j+1 j+2]) = V(:,[j+1 j+2]) / E;
% end;
% end;
% end;


%this takes care if there is a zero column.
% if ~(S(b+1,b+1)),
%     S = S(:,[1:b,(b+2):n,b+1]);
%     V = V(:,[1:b,(b+2):n,b+1]);
%     disp('hi');
% end;
% D = diag(S,1);
% end;

% Now S is diagonal. Make it nonnegative.

for j = 1:min(m,n)
if S(j,j) < 0
S(j,:) = -S(j,:);
%U(:,j) = -U(:,j);
%disp('hi');
end;
%end;

% Squeeze factors to lower right to enforce divisibility condition.

for i = 1 :  min(m,n)
for j = i+1 : min(m,n)
% Replace S(i,i), S(j,j) by their gcd and lcm respectively.
%fprintf('i is: %i, j is : %i\n',i,j)
a = S(i,i);
b = S(j,j);
if ~(a==0 && b==0),
[g,c,d] = gcd(a,b);
E = [ 1 d ; -b/g a*c/g];
F = [ c 1 ; -b*d/g a/g];
S([i j],[i j]) = E * S([i j],[i j]) * F';
%U(:,[i j]) = U(:,[i j]) / E;
%V(:,[i j]) = V(:,[i j]) / F;
%fprintf('i is: %i, j is : %i\n',i,j);
end
end;
end;

%U = round(U);
%V = round(V);
%if nargout <= 1
%U = diag(S);
end