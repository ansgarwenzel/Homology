function [B] = boundary_dihedral(degree,k,degenerate) %calculation of boundary matrices for dihedral quandles only.

if degree < 2,
    return;
end

if degenerate,
    B1 = k*((k-1)^(degree-1));
    B2 = k*((k-1)^(degree-2));
else
    B1 = k^degree;
    B2 = k^(degree-1);
end

B = sparse(B1,B2,0);

rows = cell(1,degree,1);
columns = cell(1,degree-1,1);

for i = 1:degree,
    rows{:,i,:} = 0:(k-1);
end

for i = 1:(degree-1),
    columns{:,i,:} = 0:(k-1);
end

row_names = allcomb(rows{:});
column_names = allcomb(columns{:});
[r1 r2] = size(row_names);
[c1 c2] = size(column_names);

if degenerate,
    for i = 0:(r1-1),  %delete doubles in row_names
        done=0;
        for j = 2:r2,
            if done,
            elseif ~(round(row_names(r1-i,j)-row_names(r1-i,j-1))),
                row_names(r1-i,:)=[];
                done=1;
            end
        end
    end

    for i = 0:(c1-1),  %delete doubles in column_names
        done=0;
        for j = 2:c2,
            if done,
            elseif ~(round(column_names(c1-i,j)-column_names(c1-i,j-1))),
                column_names(c1-i,:)=[];
                done=1;
            end
        end
    end
end

[r1 r2] = size(row_names);
[c1 c2] = size(column_names);

for i = 1:r1, % go through rows to calculate boundary for each row. this is what calculates the boundary!
    result_vector = zeros(1,c1);
    for j=1:r2,
        name_vector = row_names(i,:);
        b = name_vector(j);
        name_vector(j) = [];
        position = strmatch(name_vector,column_names);
        if mod(j,2),
            result_vector(position) = result_vector(position) + 1;
        else
            result_vector(position) = result_vector(position) - 1;
        end
        %use up action for quandles.
        
        if (j>1),
            for l=1:(j-1),
                name_vector(l) = up_action(name_vector(l),b,k);
            end
            for l=(j+1):(r2-1),
                name_vector(l) = down_action(name_vector(l),b,k);
            end
        end
        position = strmatch(name_vector,column_names);
        if mod(j,2),
            result_vector(position) = result_vector(position) - 1;
        else
            result_vector(position) = result_vector(position) + 1;
        end
    end
    B(i,:) = result_vector(:);
end
return;
end