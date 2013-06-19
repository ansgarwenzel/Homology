function [i,j] = rank_G(Dhat)
[m,n]=size(Dhat);
if zeros(m,n)==Dhat
    i = NaN;
    j = NaN;
    return
end
for i=1:m,
    for j=1:n,
        if Dhat(i,j) ~=0
            return;
        end
    end
end
end