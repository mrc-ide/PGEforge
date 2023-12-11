## Overview
The program takes as input a file of genotype calls for a set of samples and outputs estimates of the fraction of their genomes that are identical by descent (IBD), based on the observed genotypes and the estimated allele frequencies in the population. Optionally, a second file of genotype calls can be supplied, which will be treated as coming from a different population with different allele frequencies. For a single population, all pairwise comparisons are made between the samples (unless otherwise specified with a -b or -g flag.) For two populations, all comparisons are made between samples from the two different populations. 

Under the HMM, each variant site is assumed to be in one of two hidden states, IBD or not-IBD.  To calculate the probability of each state, estimates of the allele frequencies for every variant are required.  By default, they are calculated from the input data, but a separate file of allele frequencies can be supplied by the user (preferable if analyzing a subset of the data).

The model has two free parameters, (1) the fraction of the genome that is IBD, and (2) the number of generations during which recombination has been breaking down IBD blocks. (Note: the former is generally estimated more accurately than the latter, and is pretty robust to the latter.) The program fits for optimal values of these parameters by an iterative estimation-maximization procedure. Iterations of the fit are capped at a user-settable maximum (default = 5). To accurately determine the IBD fraction for large shared chromosome segments, only a few iterations are needed, while for smaller, older blocks of IBD, the fit may continue to improve for 15 or more iterations.

The Viterbi algorithm calculates the best single set of state assignments given data under the HMM and outputs that set. The forward-backward algorithm sums the fraction of the genome that is IBD over all possible state assignments given data under the HMM, weighting each by the probability of that set of states. If you are interested in the IBD fraction, rather than precisely which parts of the genome are IBD, this is probably the output you want (see **fract_sites_IBD** in the Output section).

It may be a good idea to thin the variant sites being supplied to hmmIBD for several reasons. First, variants with very low or zero minor allele frequency contribute little to IBD determination but do waste compute time. Second, very closely spaced sites will (be default) be automatically thinned by hmmIBD and not in an intelligent way, in an effort to avoid mutations (or sequencing errors) that extend beyond a single base pair. Third, if the density of sites varies a lot across the genome, hmmIBD's power to detect short IBD segments will be higher where site density is higher. The repository includes a Python script (thin_sites.py) that implements a simple scheme for thinning sites.

## Execution

hmmIBD is run from the command line. It requires the user to supply two options when invoked; it can also take seven optional arguments:

```
hmmIBD -i <input filename (for pop1, if using 2 pops)> -o <output filename> 
     [-I <input file, pop2>] [-f <allele frequency file (pop1)>] [-F <allele freq file (pop2)>]
     [-b <file with samples to skip>] [-m <max fit iteration>] [-n <max N generation>] [-g <file with sample pairs to use>] [-r <fixed IBD prior>]
```
Required options:
- -i: File of genotype data. See below for format.
- -o: Output file name. Two output files will be produced, with ".hmm.txt" 
      and ".hmm_fract.txt" appended to the supplied name. See below for details

Optional options:
- -f: File of allele frequencies for the sample population. Format: tab-delimited, no header, one variant per row. Line format: `<chromosome (int)> <position (bp, int)> <allele 1 freq> <all 2 freq> [<all 3 freq>] ...` The genotype and frequency files must contain exactly the same variants, in the same order. If no file is supplied, allele frequencies are calculated from the input data file.
- -I: File of genotype data from a second population; same format as for -i. (added in 2.0.0)
- -F: File of allele frequencies for the second population; same format as for -f. (added in 2.0.0)
- -m: Maximum number of fit iterations (defaults to 5).
- -b: File of sample ids to exclude from all analysis. Format: no header, one id (string) per row. (Note: b stands for "bad samples".)
- -g: File of sample pairs to analyze; all others are not processed by the HMM 	(but are still used to calculate allele frequencies). Format: no header,	tab-delimited, two sample ids (strings) per row. (Note: "g" stands for 	"good pairs".)
- - -n: Cap on the number of generations (floating point). Sets the maximum value for that parameter in the fit. This is useful if you are interested in recent IBD and are working with a population with substantial linkage disequilbrium. Specifying a small value will force the program to assume little recombination and thus a low transition rate; otherwise it will identify the small blocks of LD as ancient IBD, and will force the number of generations to be large.
- -r: Supplies a fixed value of the IBD fraction (fract_sites_IBD) used to determine IBD segments. This is useful when using hmmIBD to detect or characterize selective sweeps. Without this, IBD segments are more likely to be detected when comparing relatives, since the overall relatedness biases the probability of detecting any one segment.  

## Input file formats

Format for genotype file: tab-delimited text file, with one single nucleotide polymorphism (SNP) per line. The first two columns are the chromosome and position, followed by one sample per column. A header line, giving the sample names, is required. Genotypes are coded by number: -1 for missing data, 0 for the first allele, 1 for the second, etc. SNPs and indels (if you trust them) can thus be treated on an equal footing. The variants must be in chromosome and position order, and can have between two and eight alleles (more, if you feel like changing *max_allele* in the code).

## Output files

Two output files are produced. The file *\<filename\>.hmm.txt* contains a list of all segments (where a segment is one or more contiguous variant sites in the same state) for each sample pair, the assigned state (IBD or not-IBD) for each, and the number of variants covered by the segment; note that an assigned state of 0 means IBD, while 1 means not-IBD. These segments represent the most probable state assignments.

The file *\<filename\>.hmm_fract.txt* summarizes results for each sample pair (including some information that may be of interest only to me). If you are only interested in the fraction of the genome that is IBD between a pair of samples, look at the last column (**fract_sites_IBD**). Columns:

- **N_informative_sites**: Number of sites with data for both samples, and with at least one copy of the minor allele.
- **discordance**: Fraction of informative sites that have different alleles in the two samples.
- **log_p**: natural logarithm of the probability of the final set of state assignments and the set of observations (calculated with the Viterbi algorithm).
- **N_fit_iteration**: Number of iterations carried out in EM fitting. 
- **N_generation**: Estimated number of generations of recombination (1 of 2 free parameters in fit).
- **N_state_transition**: Number of transitions between IBD and not-IBD states across entire genome.
- **seq_shared_best_traj**: Fraction of sequence IBD based on the best state assignment, calculated as (seq in IBD segments) / (seq in IBD segments + seq in not-IBD segments). Segments in which there is a state transition between IBD and not-IBD are ignored. 
- **fract_sites_IBD**: Fraction of variant sites called IBD calculated for all possible state assignments, weighted by their probability (equal to the marginal posterior probability of the IBD state calculated with the forward-backward algorithm).
- **fract_vit_sites_IBD**: Fraction of variant sites called IBD calculated for the best state assignment (i.e. the result of the Viterbi algorithm, as in seq_shared_best_traj).


## Some variables you might want to change

The following variables control program execution in various ways, and can be easily changed in the C code prior to compilation. 

- *eps* -- assumed genotype error rate (rate of calling allele i as allele j) (default = 0.1%).
- *min_inform* -- minimum number of informative sites (those with at least one copy of the minor allele in this pair) required for processing sample pair.
- *max_discord* -- maximum fraction of informative sites with discordant genotypes between the two samples; used to skip unrelated pairs .
- *min_discord* -- minimum discordance, useful for skipping identical pairs.
- *nchrom* -- number of chromosomes in genome (14 for P. falciparum and P. vivax)
- *min_snp_sep* -- number of bp separation required between SNPs, to avoid correlated mutations spanning >1 bp (default = 10).
