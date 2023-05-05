# summary_stats

This script will calculate a set of summary statistics from multilocus sequence data. Summary statistics include nucleotide diversity (pi), number of segregating sites (SS), Watterson's theta (per locus), Watterson's theta (per site), and Tajima's D. Loci are assumed to be in nexus format (.nex or .nexus file ending) and a traits file associating individuals to species is required (see below for details). A minimum of two alleles per species are required to calculate summary statistics for any locus. Sequence data are read with Biopython (Cock et al. 2009) and summary statistics are calculated with DendroPy (Sukumaran & Holder 2010). For each species in the results file (**results_sumstats.txt**), results show the number of loci (N) used in the calculation and the mean and standard deviation of that summary statistic.

usage:  
```python
    python sumstats.py traits.file /path/to/loci/
```

***
The traits file is tab-delimited and assigns individuals to species:

ind1&nbsp;&nbsp;&nbsp;&nbsp;sp1  
ind2&nbsp;&nbsp;&nbsp;&nbsp;sp1  
ind3&nbsp;&nbsp;&nbsp;&nbsp;sp1  
ind4&nbsp;&nbsp;&nbsp;&nbsp;sp2  
ind5&nbsp;&nbsp;&nbsp;&nbsp;sp2  
ind6&nbsp;&nbsp;&nbsp;&nbsp;sp3  
ind7&nbsp;&nbsp;&nbsp;&nbsp;sp3  
ind8&nbsp;&nbsp;&nbsp;&nbsp;sp3  
...
***
References:

Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423.  

Sukumaran, J. and Mark T. Holder. 2010. DendroPy: A Python library for phylogenetic computing. Bioinformatics, 26, 1569-1571.
