# ThetaRecov
Python package to correctly estimate Tajima's D statistics and the population mutation rate from a vcf file with missing data

## What is Tajima's D
Tajima's D is the normalized difference in two emprical estimators of population mutation rate, <i>&#952;</i><sub>&#960;</sub> and <i>&#952;</i><sub>W</sub>. This is scaled to be behaving such that this would be the same as in a neutrally evolving population of constant size.
*  <i>&#952;</i><sub>&#960;</sub>: a measure of population mutation rate calculated from the average number of pairwise nucleotide differences among haplotype copies of a locus
* <i>&#952;</i><sub>W</sub>: a measure of population mutation rate calculated from the number of segragating site (S), also referred to as Watterson's <i>&#952;</i>. 


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
