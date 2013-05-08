function [B] = boundary_dihedral(degree,k,degenerate) %calculation of boundary matrices for dihedral quandles only.
if degree < 2,
    return;
end

if degenerate,
    size_matrix = [k*((k-1)^(degree-1)) k*((k-1)^(degree-2))];
    row_names = zeros(k*((k-1)^(degree-1)),degree);
    column_names = zeros(k*((k-1)^(degree-2)),degree-1).';
else
    size_matrix = [k^degree k^(degree-1)];
    row_names = zeros(k^degree,degree);
    column_names = zeros(k^(degree-1),degree-1).';
end

B = zeros(size_matrix);

repetitions_row = size_matrix(1,1)/k;
repetitions_column = size_matrix(1,2)/k;

for i=0:(k-1),
    row_names(1+i*repetitions_row:(i+1)*repetitions_row,1) = i;        
end
%the following only works for the case degree=2.
if degenerate,
    for j=2:degree,
        for i=0:(k-1);
            input_vector = 0:(k-1);
            input_vector = input_vector(input_vector~=row_names((k-1)*i+1,j-1));
            row_names(1+i*(k-1):(i+1)*(k-1),j) = input_vector;
        end
    end
else
    for j=2:degree,
        for i=0:(k-1),
            row_names(1+i*k:(i+1)*k,j) = 0:(k-1);
        end
    end
end



for i=0:(k-1),
    for j=1:(degree-1),
        column_names(1,1+i*repetitions_column:(i+1)*repetitions_column)=i;
    end
end



end