# AMP-pl-variant-selection

A collection of perl scripts for selecting _rev-env_ (RE) variants for synthesis
using a 9.5% AA frequency cutoff for gp160. A consensus _rev_ is appended to this
sequence.

IMPORTANT: This workflow will only process _rev-env-nef_ (REN) amplicons in its current
state. Processing for the _gag-pol_ (GP) amplicons may be pursued in the future.

### Setup and Usage

A snakemake framework is provided for running the scripts. The directory `data/`
should contain .fasta collections arranges like:

```
data/
   sample-1.fasta
   sample-2.fasta
    ...
```

The input .fasta files (i.e. `sample-1.fasta`, `sample-2.fasta` above) should be
single genome nucleotide sequences. Whether the input is aligned or not will not
affect the result, though you may get a warning message from BLAST regarding gap
characters. The input sequences can be uppercase or lowercase.

Finally, verify that the path to blast at the top of the `Snakefile` is correct.

```python
configfile: "config.yaml"
BLAST = "/usr/local/bin" #configure path to blast executable
DATASETS = [d for d in config for s in config[d]]
SAMPLES = [s for d in config for s in config[d]]
```

Execute by running `snakemake` via command line. The results will be output to the
`variant_analysis/` directory. The principal output files will be the .fasta file
selecting for synthesis (which end in a variant number like `*00[0-9]_s.fasta`) and
a log file ending in `*_log.txt`. There will also be a `sample_summary.csv` file for
each submitted dataset with relevant statistics for each sample.
