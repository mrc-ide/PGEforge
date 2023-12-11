## Summary sheet

| Feature           | Description   |
|-                  |-              |
|Main use           | Inference of COI, proportions, haplotypes and within-sample IBD |
|Authors            | Sha Joe Zhu, Jacob Almargo-Garcia, Gil McVean |
|Latest version     | v0.6-beta |
|Website            | https://deploid.readthedocs.io/en/latest/index.html |
|Code repository    | https://github.com/DEploid-dev/DEploid |
|Publication        | https://doi.org/10.7554/eLife.40845 (DEploidIBD), https://doi.org/10.1093/bioinformatics/btx530 (DEploid) |
|Tutorial Author    | Jason A. Hendry |
|Tutorial Date      | 11-Dec-2023 |

## Purpose
The main purpose of `DEploidIBD` is to infer complexity of infection (COI), strain proportions, strain haplotypes, and within-sample IBD from *P. falciparum* whole genome sequencing (WGS) data.

## Existing resources
- There is a [github repository](https://github.com/DEploid-dev/DEploid) containing the source code. Note that both `DEploid` and `DEploidIBD` run from the same codebase, with the later being invoked when the `-ibd` flag is used.
- There is also a [website](https://deploid.readthedocs.io/en/latest/index.html ) with additional information about installaion and running the tools.

## Citation
If you run `DEploidIBD`, please use the following citation:
```
@article{RN96,
   author = {Zhu, Sha Joe and Hendry, Jason A. and Almagro-Garcia, Jacob and Pearson, Richard D. and Amato, Roberto and Miles, Alistair and Weiss, Daniel J. and Lucas, Tim C. D. and Nguyen, Michele and Gething, Peter W. and Kwiatkowski, Dominic and McVean, Gil},
   title = {The origins and relatedness structure of mixed infections vary with local prevalence of P. falciparum malaria},
   journal = {eLife},
   volume = {8},
   pages = {e40845},
   ISSN = {2050-084X},
   DOI = {10.7554/eLife.40845},
   url = {https://doi.org/10.7554/eLife.40845},
   year = {2019},
   type = {Journal Article}
}
```

Or if you run `DEploid`, use:

```
@article{Zhu2018,
   author = {Zhu, S. J. and Almagro-Garcia, J. and McVean, G.},
   title = {Deconvolution of multiple infections in Plasmodium falciparum from high throughput sequencing data},
   journal = {Bioinformatics},
   volume = {34},
   number = {1},
   pages = {9-15},
   ISSN = {1367-4803 (Print)
1367-4803},
   DOI = {10.1093/bioinformatics/btx530},
   year = {2018},
   type = {Journal Article}
}
```