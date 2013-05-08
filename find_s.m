function [i] = find_s(S,q_minus_rho)
    for i=1:q_minus_rho
        if S(i,i)==0
            return ;
        end
    end
end