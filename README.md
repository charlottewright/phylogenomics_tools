# phylogenomics_tools
Scripts that can come in handy when analysing phylogenomic data

## test_for_base_composition_bias.py

This script performs chi-squared tests to test whether there is base composition bias across all codon positions, and at each codon position. There are two optional features:

- **RY-recoding** -  if composition bias is found at a given codon position(s), RY-recodng is performed at that position(s). The recoded alignment is then tested again for base composition bias.
- **RY-recoding and deletion** - removal of codon positions which still show bias after recoding.

The default behaviour is just to test for base composition bias, without performing recoding or deleting positions.

The default output is `output_prefix.base_composition_results.tsv`. If RY recoding is enabled, then the recoded alignment is also output `output_prefix.recoded`. If deletion is enabled, then the modified alignment is also output `output_prefix.recoded.positions.deleted`.

The last line in `output_prefix.base_composition_results.tsv` gives a handy summary of what adjustments were performed, if recoding and deletion are enabled. This is useful for deciding which modified alignment file to pick for downstream analyses.
```
test_for_base_composition_bias.py -f alignment.aln -p output_prefix # just test for bias
test_for_base_composition_bias.py -f alignment.aln -p output_prefix -r True #  test for bias and perform RY-recoding is bias is found
test_for_base_composition_bias.py -f alignment.aln -p output_prefix -r True -d True #  test for bias and perform RY-recoding is bias is found, if positions still show bias, remove them from the alignment.
```

Full usage:

```
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Path to the input file
  -p PREFIX, --prefix PREFIX
                        Prefix for the output files
  -r RECODE, --recode RECODE
                        Enable recoding of codon positions which show composition bias to RY
  -d DELETE, --delete DELETE
                        Enable removal of codon positions which show composition bias after RY recoding
```
