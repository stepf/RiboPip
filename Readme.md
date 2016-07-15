RiboPip [![Code Climate](https://codeclimate.com/github/stepf/RiboPip/badges/gpa.svg)](https://codeclimate.com/github/stepf/RiboPip) [![Dependency Status](https://gemnasium.com/badges/github.com/stepf/RiboPip.svg)](https://gemnasium.com/github.com/stepf/RiboPip) [![Inline docs](http://inch-ci.org/github/stepf/RiboPip.svg)](http://inch-ci.org/github/stepf/RiboPip)
=========
**An alignment and analysis pipeline for Ribo-seq and RNA-seq data**

RiboPip is Ruby based pipeline for processing Ribosome Profiling (Ribo-seq) and RNA sequencing (RNA-seq) datasets. `ribopip align` starts from raw sequences files and computes a splice-aware alignment to a reference database along with a read summarization (counting mapped reads per genomic feature) and data quality assessments. Read summarization data can be merged and used for differential expression estimation with `ribopip postproc`

Getting started
---------------
RiboPip wraps around all pipeline steps and ties them together for which it uses a variety of external software. It has been designed to run on every POSIX compliant UNIX system, for example, Linux, Mac OS X, and OpenBSD.

**Automatic installation**

Running `scripts/bootstrap` installs all external dependencies and sets up RiboPip correctly. It works out-of-the-box with most Linux flavours, although you might want to modify the bash script according to your needs.

**Manual installation**

1. Install all external software dependencies (see `scripts/bootstrap` for a list)
2. Manually build C extension and Ruby gem:

```bash
cd "ext/fastq_bucketize" && make # && copy bin to any directory in your PATH
cd -
rake build && gem install "./pkg/ribopip-$(rake version).gem"
```

3. Run tests:

```bash
rake spec
```

Usage
---------------
```bash
> ribopip -h
Commands:
  ribopip align
  ribopip help [COMMAND]
  ribopip postproc
```

**Example of a full pipeline run**
```bash
# Reference files:
ncrna="GRCm38.80.ncrna.fa"
genome="GRCm38.dna.primary_assembly.fa"
annotation="Mus_musculus.GRCm38.80.gtf"
igv="GRCm38.igv.genome"
# Arguments for ribopip align:
aln_args="-n ${ncrna}$ -g ${genome} -a ${annotation} --igv_ref ${igv}"

# Run upstream pipeline for datasets, each comprising of
#  - footprint and mrna data   : fp, mrna
#  - 2 experimental conditions : treated, control
#  - 2 replicates              : rep1, rep2
ribopip align $aln_args -r fp_treated_rep1.fastq
ribopip align $aln_args -r fp_treated_rep2.fastq
ribopip align $aln_args -r fp_control_rep1.fastq
ribopip align $aln_args -r fp_control_rep2.fastq
ribopip align $aln_args -r mrna_treated_rep1.fastq
ribopip align $aln_args -r mrna_treated_rep2.fastq
ribopip align $aln_args -r mrna_control_rep1.fastq
ribopip align $aln_args -r mrna_control_rep2.fastq

# Run downstream pipeline to analyze feature counts
ribopip postproc -a $annotation -o . \
  --fp_1 fp_treated \
  fp_treated_rep1.vs_genome.uni.ft.dist.txt \
  fp_treated_rep2.vs_genome.uni.ft.dist.txt \
  --fp_control fp_control \
  fp_control_rep1.vs_genome.uni.ft.dist.txt \
  fp_control_rep2.vs_genome.uni.ft.dist.txt \
  --mrna_1 mrna_treated  \
  mrna_treated_rep1.vs_genome.uni.ft.dist.txt \
  mrna_treated_rep2.vs_genome.uni.ft.dist.txt \
  --mrna_control mrna_control  \
  mrna_control_rep1.vs_genome.uni.ft.dist.txt \
  mrna_control_rep2.vs_genome.uni.ft.dist.txt
```

## Copyright
Copyright (c) 2016 Stefan Dang. See LICENSE for details.
