# SIMD-package
Statistical Inferences with MeDIP-seq Data (SIMD) to infer the methylation 
level for each CpG site.

# Introduction
DNA methylation is an essential epigenetic modification involved in regulating 
the expression of mammalian genomes. A variety of experimental approaches to 
generate genome-wide or whole-genome DNA methylation data have emerged in 
recent years. Methylated DNA immunoprecipitation followed by sequencing 
(MeDIP-seq) is one of the major tools used in whole-genome epigenetic studies. 
However, analyzing this data in terms of accuracy, sensitivity, and speed still 
remains an important challenge. Existing methods, such as BATMAN and MEDIPS, 
analyze MeDIP-seq data by dividing the whole genome into equal length windows 
and assume that each CpG of the same window has the same methylation level. 
More precise work is necessary to estimate the methylation level of each CpG 
site in the whole genome. In this paper, we propose a Statistical Inferences 
with MeDIP-seq Data (SIMD) to infer the methylation level for each CpG site. 
This package provides a inferential analysis method for detecting differentially 
expressed CpG sites in MeDIP-seq data. It uses statistical framework and EM 
algorithm, to identify differentially expressed CpG sites.

#User's Guide
Please refer to the "SIMD.pdf" vignetee for detailed function 
instructions.
[![Travis-CI Build Status](https://travis-ci.org/FocusPaka/SIMD.svg?branch=master)](https://travis-ci.org/FocusPaka/SIMD)