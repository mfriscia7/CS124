function [ variant_count ] = findfreq( Reads, StrainSet, StrainSetME, location)

%%Find a way to map reads to strain
SNP_num = size(location,2);
variant_count = zeros(SNP_num,1);
read_len = size(Reads,2);
strain_len = size(StrainSet,2);


%%Make a matrix of SNP strings with length of read on either side
SNP_check = zeros(SNP_num,(read_len*2)-1);

%%for each SNP
for i=1:SNP_num
   for j=1:(read_len*2)-1
      if (location(1,i) - read_len + j > 0 && location(1,i) - read_len + j <= strain_len)
        SNP_check(i,j) = StrainSet(1, location(1,i) - read_len + j);
      end
   end
   SNP_check(i,read_len) = location(2,i);
end



%%this is the trivial mapping but it works better because it's not tested
%%against the whole genome, just the areas near the SNPs

%%check if SNPs are found in each Read
for i=1:size(Reads,1)
  %%check each SNP
  for j=1:SNP_num
     %%check each position
     for k=1:read_len
         if (Reads(i,:) == SNP_check(j,k:k+read_len-1))
            variant_count(j,1) = variant_count(j,1) + 1; 
         end
     end
  end
end

%%turn variant count into percentages
variant_count = variant_count / sum(variant_count);



%%Least-Squares
answer = lsqnonneg(StrainSetME,variant_count);


end

