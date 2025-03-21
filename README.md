# ThetaRecov
Python package to correctly estimate Tajima's D statistics and the population mutation rate from a vcf file with missing data

## What is Tajima's D
Tajima's D (Tajima 1989) is the normalized difference in two emprical estimators of population mutation rate, <i>&#952;</i><sub>&#960;</sub> and <i>&#952;</i><sub>W</sub>. This is scaled to be behaving such that this would be the same as in a neutrally evolving population of constant size. The equation is as follows:  
```math
D = \frac{\hat{\theta}_{\pi} - \hat{\theta}_{W}}{\sqrt{\mathrm{Var}(\hat{\theta}_{\pi} - \hat{\theta}_{W})}}
```

*  <i>&#952;</i><sub>&#960;</sub>: a measure of population mutation rate calculated from the average number of pairwise nucleotide differences among haplotype copies of a locus
* <i>&#952;</i><sub>W</sub>: a measure of population mutation rate calculated from the number of segragating site (<i>S</i>), also referred to as Watterson's <i>&#952;</i>. 


## Installation
To install:
```bash
pip install git+https://github.com/akihirao/ThetaRecov.git
```

To update:
```bash
pip install --upgrade git+https://github.com/akihirao/ThetaRecov.git
```

## Working examples
```bash
# Computation of theta and Tajima's D across overall genome
tajimaD_overall --input_vcf test.vcf.gz --output_csv test_TajimaD_overall.csv

# Computation of theta and Tajima's D in sliding windows
tajimaD_windows --input_vcf test.vcf.gz --output_csv test_TajimaD_windows.csv --windows_size 1000
```

## Reference
Tajima F (1989) Statistical Method for Testing the Neutral Mutation Hypothesis by DNA Polymorphism. Genetics, 123: 585–595. https://doi.org/10.1093/genetics/123.3.585
