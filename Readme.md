RiboPip [![Code Climate](https://codeclimate.com/github/stepf/RiboPip/badges/gpa.svg)](https://codeclimate.com/github/stepf/RiboPip) [![Dependency Status](https://gemnasium.com/badges/github.com/stepf/RiboPip.svg)](https://gemnasium.com/github.com/stepf/RiboPip) [![Inline docs](http://inch-ci.org/github/stepf/RiboPip.svg)](http://inch-ci.org/github/stepf/RiboPip)
=========
**An alignment and analysis pipeline for Ribo-seq and RNA-seq data**

RiboPip is Ruby based pipeline for processing Ribosome Profiling (Ribo-seq) and RNA sequencing (RNA-seq) datasets. `ribopip align` starts from raw sequences files and computes a splice-aware alignment to a reference database along with a read summarization (counting mapped reads per genomic feature) and data quality assessments. Read summarization data can be merged and used for differential expression estimation with `ribopip postproc`

Getting started
---------------
RiboPip wraps around all pipeline steps and ties them together for which it depends a variety of external software. It has been designed to run on every POSIX compliant UNIX system, for example, Linux, Mac OS X, and OpenBSD.

**Automatic installation**

Running `scripts/bootstrap` installs all external dependencies and sets up RiboPip correctly. It works out-of-the-box with most Linux flavours, although you might want to modify the bash script according to your needs.

**Manual installation**

* Install all external software dependencies (see `scripts/bootstrap` for a list)
* Manually build C extension and Ruby gem:

```bash
cd "ext/fastq-bucketize-0.1" && make # && copy bin to any directory in your PATH
cd -
rake build && gem install "./pkg/ribopip-$(rake version).gem"
```

* Run tests:

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

Pipeline feature overview
---------------

### Pre-Processing: Mapping the sequencing results
The pre-processing pipeline comprises of four consecutive steps. The first step filters and prepares the raw sequence files. The second step removes unwanted RNAs. The third step computes a splice-aware alignment to a given reference database. The fourth step extracts data subsets with desired properties for further analysis.

1. **Data preparation**
  * **Linker clipping**: Deep-sequencing technologies require specific linker sequences to be ligated to the 3’ fragment ends, which can introduce an analytical bias. Users can provide a linker sequence and choose between an perfect-match or error-tolarant clipping approach.
  * **Nucleotide trimming**: During reverse transcription an untemplated nucleotide is frequently added to the 5’ end of each read, which can be trimmed off.
  * **Length selection**: Very short reads can introduce ambiguities and thus can be removed.
2. **Removal of unwanted RNAs**: Nuclease footprinting routinely leaves rRNA intact, which can comprise a large fraction of the reads and introduce a bias. All unwanted RNAs can be removed by mapping them contiguously to a corresponding sequence database (e.g. all ncRNAs) and discarding all successfully mapped reads.
3. **Alignment to the genome reference**: The remaining reads can be aligned to a genome reference using a splice-aware aligner.
4. **Read extraction**: For more specific analysis, subsets of reads can be extracted based on criteria like the number of hits (unique, multiple) or number of mismatches.

#### Metrics
During pre-processing a number of metrics is computed:

* **Read counts**: Counting the total number of alignments contained in a file provides a simple, yet important metric about a pre-processing step, e.g. how many reads were successfully mapping to the ncRNA database and sorted out consequently.
* **Counting reads per genomic feature**: The alignment to a genome reference during pre-processing results in reported coordinates for each aligned read. These coordinates can be further matched to known annotated genomic features. This enables to test for Differential Expression or Translational Efficiency.
* **Quality metrics**: Quality checking reads can help to spot potential issues before the alignment is performed. For example a dropping base quality towards the 3’ ends of all reads can be typical for NGS data.
* **Density tracks**: Much information in a ribosome profiling experiment is contained in the footprint density. Visual inspection using a genome browser programs like IGV enable investigators to not only quality control, but furthermore to inspect potential single nucleotide variants (SNVs) or alternative splicing sites (Sashimi Plots) at all sites of interest.

### Post-Processing: Statistical analyses

* **Expression normalization**: By assigning reads to annotated genomic features, it is possible to estimate the expression of that feature. As read counts arising from a transcript are proportional to the transcript length and sampling depth, it is important to normalize these expression estimates to ensure comparability across different samples and features.
* **Differential expression**: Information across experimental replicates is used to estimate normalization factors and dispersions to then test each genomic feature for differential expression.
* **Translational efficiency**: Information about ribosome occupancy and RNA abundance can be exploited to determine the translational efficiency every gene.

Copyright
---------------
Copyright (c) 2016 Stefan Dang. See LICENSE for details.
