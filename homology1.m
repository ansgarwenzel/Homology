% function [Delta,k,s_minus_k] = homology1(degree,l,degenerate)
function [] = homology1(degree,l,degenerate)

F=boundary_dihedral(degree+1,l,degenerate);
G=boundary_dihedral(degree,l,degenerate);
[~,q1]=size(F);
[q,~]=size(G);
if q~=q1
    Disp('help, something is wrong')
    Delta=NaN;
    k=NaN;
    s_minus_k=NaN;
    return;
end


%calculate X and Y
[~,Y] = Hermite_normal_form(G);
[~,X] = Hermite_normal_form((G*Y).');

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

S=smith(N);

s=find_s(S,q_minus_rho);

Delta=sparse(s,s,0);
s_minus_k=0;
k=0;
for i=1:s,
    Delta(i,i)=S(i,i);
    if Delta(i,i)==1,
        s_minus_k=i;
    else
        k=k+1;
    end
end
Delta = full(Delta);
%display_vector =  repmat('%i ',1,s);
other_terms = diag(Delta);
other_terms = other_terms(((s_minus_k+1):end));
%fprintf('the Matrix Delta is: \n\n');
%for i=1:s,
%    fprintf(strcat(num2str(display_vector) , ' \n' ),full(Delta(i,1:s)));
%end
fprintf('There are %i ones and %i other terms\n',s_minus_k,k);
fprintf('The other terms are:\n');
for i=1:k,
    fprintf('%i\n',other_terms(i));
end
%fprintf('\n');
%fprintf('And then H_%i = Z^%i + %i other terms.\n',degree,s_minus_k,k);

%return;
end