# DistKnKs
Calculates Kn/Ks distances from two CDS sequences.

## Description
This program will read two aligned CDS sequences (see [CodonAlign](https://github.com/santiagosnchez/CodonAlign)) and measure synonymous (Ks) and nonsynonymous (Kn) distances applying K80 and JC correction for hidden mutations.

A nice feature is that you can specify any codon table by using the `--codon_table` argument, which will be retrieved from the [NCBI translation table website](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1).
