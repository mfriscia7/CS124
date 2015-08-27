%makeSample(strain_num, strain_len, SNP_num, SampleSize)
[Sample, StrainSet, StrainSetME, location, freq] = makeSample(10,100,10,10000);

%SampleReads(Sample, read_num, read_len)
Reads = SampleReads(Sample,1000,30);

%baseline = baseline(Reads, StrainSet);
variant_count = findfreq(Reads, StrainSet, StrainSetME, location);