# get_sequence_from_vcf.py

This script extracts the most probable sequence of all genotypes of a specified vcf file list. Only the major allele is considered, there is no phasing ! This tool is reliable for haploid genomes, not for polyploid ones. 

## Install

`wget https://github.com/nlapalu/misc/archive/master.zip`

`unzip master.zip`

`cd misc-master`

python module [pysam](https://github.com/pysam-developers/pysam) must be installed.

## Usage and options

### Usage:

`python get_sequence_from_vcf.py "chr_1" 1000 2000 vcf.list -v 2`

### Arguments:

| Argument | Description |
| --------- | ----------- |
| contig | Sequence name | 
| start | Start coordinate, 1-based indexing system |
| end | End coordinate, 1-based indexing system |
| listFH | file with list of bgzip vcf files to analyze, one file per line |

### Options:

| Option | Description |
| ------ | ----------- |
| `-f, --freq` |  minimal frequence to consider the allele, filtered under this value [default=0.0] |
| `-n, --nbAlleles` | maximal number of alleles allowed per sample, filtered upper this value [default=2] |
| `-mi, --minDP` | minimal depth of coverage to consider the allele, filtered under this value, [default=0] |
| `-ma, --maxDP` | maximal depth of coverage to consider the allele, filtered upper this value, [default=0] |
| `-sf, --seqFasta` | export sequence in fasta file |
| `-af, --algmtFasta` | export alignment in fasta file |
| `-v, --verbosity` | increase output verbosity 1=error, 2=info, 3=debug |
| `--version` | tool version |
| `-h, --help` | help message |


## Inputs

Inputs are vcf files compressed in *Blocked GNU Zip Format* with bgzip (available with tabix). Then each *\*.vcf.gz* file to analyze must be written in a file (one file per line) which is supplied to the program.

## Outputs

The output is separated in two sections: the first is a multi-fasta with all genotypes, the second is a fasta-like alignment format including the reference sequence. By default outputs are written to stdout, but you can redirect them to files with *-sf* and *-af* options.

In case of the vcf file does not contain any information on a requested position (example no mapping = no VCF entry), the program returns a "N". In the same way, if a position does not meet the required criteria (depth, frequency), the program also returns a "N".
Indels are left aligned as in alignment and vcf files, so to obtain the most reliable fasta alignment we recommend to perform a re-alignment with a dedicated tool.


#### Example 1: one deletion in two genomes (G1,G2) compare to the reference

`python get_sequence_from_vcf.py "chr_1" 818 831 vcf.list -v 2`

```
### sequence in fasta format ###

>G1:chr_1:818:831
TCTGA

>G2:chr_1:818:831
TNTGA

>G3:chr_1:818:831
TCTTTGTATGTGGA

>G4:chr_1:818:831
TCTTTGTATGTGGA


### alignment in fasta format ###

>reference:chr_1:818:831
TCTTTGTATGTGGA

>G1:chr_1:818:831
TCTG---------A

>G2:chr_1:818:831
TNTG---------A

>G3:chr_1:818:831
TCTTTGTATGTGGA

>G4:chr_1:818:831
TCTTTGTATGTGGA
```

#### Example 2: one insertion in two genomes compare (G1,G2) to the reference 

`python get_sequence_from_vcf.py "chr_1" 2829 2832 vcf.list -v 2`

```
### sequence in fasta format ###

>G1:chr_1:2829:2832
CTAAGGAAATAT

>G2:chr_1:2829:2832
CTAAGGAAATAT

>G3:chr_1:2829:2832
CTAT

>G4:chr_1:2829:2832
CTAT


### alignment in fasta format ###

>reference:chr_1:2829:2832
CT--------AT

>G1:chr_1:2829:2832
CTAAGGAAATAT

>G2:chr_1:2829:2832
CTAAGGAAATAT

>G3:chr_1:2829:2832
CT--------AT

>G4:chr_1:2829:2832
CT--------AT
```
