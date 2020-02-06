# AMP-pl-variant-selection

A collection of perl scripts for selecting _rev-env_ (RE) variants for synthesis
using a 9.5% AA frequency cutoff for gp160. A consensus _rev_ is appended to this
sequence.

### Setup and Usage

A snakemake framework is provided for running the scripts. The directory `postproc/`
should contain the following subdirectory structure:

```
postproc/
  dataset-1/
    sample-1/
      sample-1.fasta
    sample-2/
      sample-2.fasta
  dataset-2/
    sample-1/
      sample-1.fasta
    ...
```

If you're running the [porpid-postproc](https://gitlab.com/MurrellLab/porpid-postproc) workflow you can also just copy over the `postproc/`
directory.

Next, specify what samples you want to run in the `config.yaml` file, i.e

```yaml
dataset-1:
  - sample-1
  - sample-2
dataset-2:
  - sample-1
```

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
