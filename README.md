# ThetaRecov
Python package for correct computation of <i>&#952;</i> and Tajima's D in population genetics from a vcf file with missing data

## What is Tajima's D
Tajima's D is the normalized difference of two population genetics parameters, the average number of pairwise differences, <i>&#952;</i><sub>&#960;</sub>, and Watterson's theta, <i>&#952;</i><sub>W</sub>. This is scaled to be behaving such that this would be the same as in a neutrally evolving population of constant size.

## Requirements
* Python packages
	* cyvcf2: a fast VCF parser https://brentp.github.io/cyvcf2
	* numpy: a fundamental package for scientific computing https://github.com/numpy/numpy
	* pandas: python data analysis library https://pandas.pydata.org


## Installation
To install:
```bash
pip install git+https://github.com/akihirao/ThetaRecov.git
```

To update:
```bash
pip install --upgrade git+https://github.com/akihirao/ThetaRecov.git
```

## An example of execution
```bash
# Computation of theta and Tajima's D across overall genome
tajimaD_overall test.vcf.gz test_TajimaD_overall.csv

# Computation of theta and Tajima's D in sliding windows
tajimaD_windows test.vcf.gz 1000 test_TajimaD_windows.csv
```
