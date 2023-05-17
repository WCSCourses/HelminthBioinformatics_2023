# Genetic variation module - cheat sheet

#### 2.2. Questions:
- How does the "all sample" report compared to the "Australian-only sample" report?
- From the "General statistics" section, can we see any sample groups that look different and that might be problematic for our analyses?
- Are all of the read lengths the same? How can you tell?
- Can you think of any other quality control checks on the raw data that you might want to do before analysing your data? Hint: think about where these worms are collected from. 


#### 3.1. Questions:
- Looking at the WormBase ParaSite website for *Haemonchus contortus* - how big is the genome? how many genes are present?
- There is a second *Haemonchus contortus* genome resource also present - how do the two genomes compare?


#### 5.1. Question: 
- can you think of a reason why variant calling in a haploid sample is less complicated than compared to a diploid (or polyploid) sample?


#### 8.1. Questions:
- how many variants were kept after filtering?
- it is possible that not all samples will contain all of the SNPs, ie. there is some degree of "missingness". Can you find a flag in the vcftools manual to test this? Are there samples with a lot of missing data?


#### 9.2. Questions:
- what proportion of SNPs are in coding regions vs non-coding regions? Why would this happen?
- what proporiton of variants are a "synonymous_variant" and what proportion are "missense_variants"? What effect to these variants have on the coding sequence?
- find a gene with a missense variant - what is the amino acid change, and is it likely to have an effect on the protein? (use the following table to help you: [table](https://en.wikipedia.org/wiki/File:ProteinogenicAminoAcids.svg)) 
- can you think of other ways to determine if this variant might impact the function of this protein?


#### 11.1. Questions: 
- How do these plots compare? 
- What is the relative contribution of variance in the PC3/PC4 plot compared to the PC1/PC2 plot?


#### 11.2. Questions:
- Looking at the ellipses specifically, can you see any countries that have a different distribution than the others, and describe this difference?


#### 12.1. Questions:
- How do the samples cluster on the tree?
- Are the PCAs easier or harder to interpret than the tree? Why?