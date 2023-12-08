# Whole genome sequencing (WGS) data for *P. falciparum*

## Overview
This folder holds several WGS datasets for *P. falciparum* designed to be used as part of tutorials.

## Pf3k
All of the data within the subfolder `wgs/pf3k` was derived from the [Pf3k Project](https://www.malariagen.net/parasite/pf3k). Currently, there are three VCF files, with corresponding CSVs containing metadata, for samples from:
- Democratic Republic of the Congo ($n=113$)
- Vietnam ($n=97$)
- *In vitro* mixures of laboratory strains ($n=25$)

Each VCF contains 247,496 high-quality (VQSLOD>6) biallelic SNPs across all fourteen somatic chromosomes. The VCFs are sorted and an index file is provided. The Fws statistics provided in the metadata CSVs were collected from the [Pf7 data set](https://www.malariagen.net/sites/default/files/Pf7_fws.txt), which contains the Pf3k samples. These were not calculated for the *in vitro* lab mixtures.

