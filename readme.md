# UniParental Disomy triO investiGator (UPDOG)

A Python tool for detecting and analyzing UniParental Disomy (UPD) events in Next Generation Sequencing (NGS) trio data.

## Overview

UniParental Disomy (UPD) occurs when an individual inherits both copies of a chromosome from a single parent and none from the other. UPD manifests in two forms:
- Heterodisomy: Inheritance of non-identical chromosomes from one parent
- Isodisomy: Inheritance of a duplicate chromosome from one parent

UPDOG analyzes VCF files containing trio data (child and both parents) to identify and visualize potential UPD events through:
1. Calculating chromosome segment metrics
2. Detecting statistically significant Mendelian errors
3. Generating visualizations and detailed reports

## Installation

From PyPI:
```bash
pip install pyupdog
```

From GitHub (a branch, tag or commit can be given after the `@`):
```bash
pip install git+https://github.com/AWGL/UPDOG.git@main
```

From a local clone (use `-e` for an editable install when developing):
```bash
git clone https://github.com/AWGL/UPDOG.git
cd UPDOG
pip install .
```

All of these install the `upd` package and put `UPDog.py` on your `PATH`. Requires Python >= 3.11.

## Usage

Basic command:
```bash
UPDog.py --vcf <input.vcf.gz> --proband_id <sample_id> --ped <family.ped> --output <output_prefix>
```

Full options:
```bash
UPDog.py \
  --vcf <input.vcf.gz>           # Input VCF (must be bgzipped and tabixed)
  --proband_id <sample_id>       # Sample ID of the proband in VCF
  --ped <family.ped>             # PED file describing family relationships
  --output <output_prefix>       # Output name prefix
  --min_dp <int>                 # Minimum genotype depth (default: 20)
  --min_gq <int>                 # Minimum genotype quality (default: 20)
  --min_qual <int>               # Minimum QUAL value (default: 90)
  --block_size <int>             # Block size for UPD calculation (default: 1000000)
  --min_variants_per_block <int> # Minimum variants per block (default: 100)
  --p_value <float>              # Max P-value for significance (default: 0.001)
  --min_blocks <int>             # Minimum contiguous blocks for calls (default: 5)
  --min_proportion <float>       # Minimum UPD variant proportion (default: 0.01)
  --chromosome <str>             # Restrict to single chromosome (testing only)
  --prop_plot <bool>             # Plot variant proportions aka BAF plot (default: False)
  --wes <bool>                   # BAF plot not downsampled - better for WES data (default: False)
```

## Analysis Method

1. **Segmentation**: Examines each chromosome in fixed-size segments (default 1MB)
2. **Metric Calculation**: Computes proportion of Mendelian errors per segment
3. **Baseline Determination**: Calculates expected Mendelian error rate across all chromosomes
4. **Statistical Analysis**: Identifies segments with significant increases in Mendelian errors
5. **Merging**: Combines adjacent significant segments
6. **Filtering**: Applies quality filters to generate final UPD calls

## Output Files

- `*_UPD_calls.csv`: Tab-separated file containing significant UPD events
- `*_raw_data.csv`: Raw metrics for all analyzed segments
- `*_chr*_UPD.png`: UPD metric plots per chromosome
- `*_chr*_baf.png`: B-allele frequency plots (if --prop_plot enabled)

## Unit tests

Unit tests covering utility functions can be run as follows:
```
python upd/tests.py
```

## References

1. [UniParental Disomy - Wikipedia](https://en.wikipedia.org/wiki/Uniparental_disomy)
2. [UPDio: A Tool for Detection of Iso- and Heterodisomy in Parent-Child Trios Using SNP Microarrays](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3975066/)
