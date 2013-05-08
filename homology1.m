function [Delta,k,s_minus_k] = homology1()

F=[1,-1,0,0,1,0;1,-1,-1,0,0,0;-1,1,1,0,0,0;-1,1,0,0,-1,0;0,0,1,-1,0,1;-1,0,1,-1,0,0;0,0,-1,1,0,-1;1,0,-1,1,0,0;0,-1,0,0,1,-1;0,0,0,1,1,-1;0,0,0,-1,-1,1;0,1,0,0,-1,1];
G=[-1,0,1;-1,1,0;0,-1,1;1,-1,0;0,1,-1;1,0,-1];
[p,q1]=size(F);
[q,r]=size(G);
if q~=q1
    Disp('help, something is wrong')
    Delta=NaN;
    k=NaN;
    s_minus_k=NaN;
    return;
end


%calculate X and Y
[H Y] = Hermite_normal_form(G);
[H_2 X] = Hermite_normal_form((G*Y).');

X = X.';
X = inside_out(X);
Y = inside_out(Y);
Dhat = X*G*Y;

q_minus_rho = rank_G(Dhat)-1;

rho = q-q_minus_rho;

Identity_matrix = eye(q);
for i = (q_minus_rho+1):q
    Identity_matrix(:,i)=0;
    Identity_matrix(i,:)=0;
end

Z = Identity_matrix*X;
Z = Z(1:q_minus_rho,:);

N=round(F/Z);

[U S T]=smith(N);

s=find_s(S,q_minus_rho);

Delta=zeros(s);
s_minus_k=0;
k=0;
for i=1:s
    Delta(i,i)=S(i,i);
    if Delta(i,i)==1
        s_minus_k=i;
    else
        k=k+1;
    end
end

return;
end