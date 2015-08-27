function [ Sample, StrainSet, StrainSetME, location, freq ] = makeSample( strain_num, strain_len, SNP_num, SampleSize )

SNP_freq = 4;
initial_strain = zeros(1,strain_len);
StrainSet = zeros(strain_num+1,strain_len);
%%these two are used to make sure all SNPs are found in at least one strain
%%of the sample
SNP_added = zeros(1, SNP_num);
SNP_temp = zeros(1, SNP_num);

%%make initial strain
for i=1:strain_len
    initial_strain(1,i) = randi(4);
end

new_strain = initial_strain;



%%create SNP locations
location = zeros(2,SNP_num);
pool = randperm(strain_len);

for i=1:SNP_num 
   %%takes first few of the randperm
   location(1,i) = pool(1,i);
   
   %%determines base change
   base = randi(4);
   while base == initial_strain(1,location(1,i)); 
    base = randi(4);
   end
   location(2,i) = base;
end

SNP_count = zeros(strain_num,1);
StrainSet(1,:) = initial_strain;

while prod(SNP_added) == 0
    SNP_added = zeros(1, SNP_num);
for i=1:strain_num
    %%create a way to make SNPs in strains
    while ismember(new_strain, StrainSet, 'rows')
        SNP_temp = zeros(1, SNP_num);
      for j=1:SNP_num
        if randi(SNP_freq) == 1
            x = randi(SNP_num);
            if new_strain(1,location(1,x)) ~= location(2,x)
               new_strain(1,location(1,x)) = location(2,x);
               SNP_temp(1,x) = 1;
            end
        end
      end
    end
    
    StrainSet(i+1,:) = new_strain;
    new_strain = initial_strain;
    SNP_added = SNP_added + SNP_temp;
end
end

%%keep track of SNP number
for i=1:strain_num
    for j=1:strain_len
        if StrainSet(1,j) ~= StrainSet(i+1,j)
           SNP_count(i) = SNP_count(i) + 1; 
        end
    end
end


%%set frequencies we want to eventually find in the sample
%%summing all of the frequencies equals 1
freq = round(diff([0;sort(rand(strain_num-1,1));1])*100);

%%makes sure to get integers summing to 100
while sum(freq) ~= 100
    freq = round(diff([0;sort(rand(strain_num-1,1));1])*100);
end
freq = freq/100;
freq = sort(freq,'descend');


StrainSet_clone = StrainSet;

%%need to sort StrainSet
for i=1:strain_num
   [y,slot] = max(SNP_count);
   
   StrainSet(i+1,:) = StrainSet_clone(slot+1,:);
   SNP_count(slot) = 0;
end

%%create Sample
Sample = zeros(SampleSize,strain_len);
counter = 1;

for i=1:size(freq,1)
    for j=1:round(freq(i,1)*SampleSize)
       Sample(counter,:) = StrainSet(i+1,:);
       counter = counter + 1;
    end
end

%%make a strain set with 1's, 0's
StrainSetME = zeros(strain_num,strain_len);
for i=2:strain_num+1
   for j=1:strain_len
      if StrainSet(1,j) ~= StrainSet(i,j)
         StrainSetME(i-1,j) = 1;
      end
   end
end

end


