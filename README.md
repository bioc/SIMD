# SIMD-package
This package use Statistical Inferences with MeDIP-seq Data (SIMD) to infer the 
methylation level for each CpG site. The methods on this package are described 
in the article 'Methylation-level Inferences and Detection of Differential 
Methylation with Medip-seq Data'[1]. 

# Introduction
DNA methylation is an essential epigenetic modification involved in regulating 
the expression of mammalian genomes. A variety of experimental approaches to 
generate genome-wide or whole-genome DNA methylation data have emerged in 
recent years. Methylated DNA immunoprecipitation followed by sequencing 
(MeDIP-seq) is one of the major tools used in whole-genome epigenetic studies. 
However, analyzing this data in terms of accuracy, sensitivity, and speed still 
remains an important challenge. Existing methods, such as BATMAN[2] and MEDIPS[3], 
analyze MeDIP-seq data by dividing the whole genome into equal length windows 
and assume that each CpG of the same window has the same methylation level. 
More precise work is necessary to estimate the methylation level of each CpG 
site in the whole genome. In this paper, we propose a Statistical Inferences 
with MeDIP-seq Data (SIMD) to infer the methylation level for each CpG site. 
This package provides a inferential analysis method for detecting differentially 
expressed CpG sites in MeDIP-seq data. It uses statistical framework and EM 
algorithm, to identify differentially expressed CpG sites.

# User's Guide
Please refer to the "SIMD.pdf" vignetee for detailed function 
instructions.

# Reference
1. Zhou Y, Zhu JD, Zhao MT, Zhang Bx, Jiang CF, and Yang XY. "Methylation-level 
Inferences and Detection of Differential Methylation with Medip-seq Data". 
Plos one. Accepted.
2. Down TA, Rakyan VK, Turner DJ, Flicek P, Li H, Kulesha E, et al. "A Bayesian
deconvolution strategy for immunoprecipitation-based DNA methylome analysis".
Nature Biotechnology. 2008; 26: 779-785.
3. Chavez L, Jozefczuk J, Grimm C, Dietrich J, Timmermann B, Lehrach H, et al.
"Computational analysis of genome-wide DNA methylation during the
differentiation of human embryonic stem cells along the endodermal lineage".
Genome Research. 2010; 20: 1441-1450.
[![Travis-CI Build Status](https://travis-ci.org/FocusPaka/SIMD.svg?branch=master)](https://travis-ci.org/FocusPaka/SIMD)
