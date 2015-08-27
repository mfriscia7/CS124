function [ Reads ] = SampleReads( Sample, read_num, read_len )

Reads = zeros(read_num, read_len);
x = zeros(1, read_len);

for i=1:read_num
    row = randi(size(Sample,1));
    column = randi(size(Sample,2) - read_len);
    
    for j=1:read_len
        x(1,j) = Sample(row, column + j);
    end
    Reads(i,:) = x;

end

