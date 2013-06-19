function [Y] = inside_out(X)
[n n] = size(X);
Y = X;
for i = 1:n,
    for j = 1:n,
        Y(n-j+1,n-i+1) = X(j,i);
    end
end
return;
end