function [ strain_count ] = baseline(Reads, StrainSet)

strain_count = zeros(size(StrainSet,1)-1,1);
count = 0;
found = 0;
read_len = size(Reads,2);
strain_len = size(StrainSet,2);

%%each Read
for i=1:size(Reads,1)
    %%test against each Strain
    for j=2:size(StrainSet,1)
        %%test against each position
        while found == 0 && count+read_len <= strain_len
            count = count + 1;
            if Reads(i,:) == StrainSet(j, count:count+read_len-1)
                found = 1;
                strain_count(j-1,1) = strain_count(j-1,1) + 1;
            end
        end
        count = 0;
        found = 0;
    end
end

strain_count = strain_count / sum(strain_count);

end

