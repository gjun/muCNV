# muCNV

Multi-sample SV genotyper for large-scale WGS data

Version 1.0.0
Last edited: July 9, 2020, (c) Goo Jun

muCNV uses multiple steps for multi-sample SV genotyping, to handle large number of samples and to enable efficient parallelization:

0. Generate binary interval file from input VCF 
1. Generate *pileups* from CRAM/BAM files using candidate list of SVs
2. Merge individual pileups to generate 'merged' pileups
3. Multi-sample genotyping from merged pileups

* The 'merging' step is optional, but is recommended for ~ >500 samples. 
* By default, merging will create pileups separated by chromosomes. 
* When genotyping merged pileups, do not add 'chr' postfixes in your 'list' file

## Installation
Required: g++ compiler environment and HTSlib (https://github.com/samtools/htslib).
If HTSlib is not installed system-wide, please manually add -I(HTSlib include directory) and -L(HISlib library directory) to Makefile. 
```
$ make
```

To test, 
```
$ bin/muCNV  
```
You can copy muCNV binary to one of your $PATH locations (e.g. /usr/local/bin).

## Basis usage:
```
$ mucnv [command] [options]
$ mucnv [command] --help
```    

# Use-case Example
1. Input files 
- GRCh38 CRAM files 
```
sample1.cram 
sample2.cram
sample3.cram
...
sample10.cram
```
- VCF file containing candidate SV variants
```
candidate_sv.vcf
```

2. Generate 'interval' file from VCF
```
$ muCNV vcf2int -v candidate_sv.vcf -i candidate_sv.interavls 
```

3. Generate pileup from each CRAM file
Assuming muCNV source directory is at /home/user1/muCNV and pileup files are being created under $(pwd)/pileup
```
$ mkdir pileup
$ muCNV pileup -s sample1 -V candiate_sv.intervals -f /home/user1/muCNV/resource/GRCh38.v3.gc -o pileup/sample1
$ muCNV pileup -s sample2 -V candiate_sv.intervals -f /home/user1/muCNV/resource/GRCh38.v3.gc -o pileup/sample2
...
$ muCNV pileup -s sample10 -V candiate_sv.intervals -f /home/user1/muCNV/resource/GRCh38.v3.gc -o pileup/sample10
```
You can simplify/parallelize pileup step by using GNU Parallel or GNU Make.

4. Genotype candidate SVs from pileup data
Make a text file containing the base names of each pileup file, say, pileup.list
```
$ cat pileup.list
pileup/sample1
pileup/sample2
...
pileup/sample10
```
Then run genotyping step to create the output file, genotyped.vcf:
```
$ muCNV genotype -V candidate_sv.intervals -i pileup.list -o genotyped.vcf -f /home/user1/muCNV/resource/GRCh38.v3.gc
```

# Commands in more details
## Generating variant list from input VCF
```
$ muCNV vcf2int -v <VCF> -i <interval file> -p

-i, --interval
  Name of the interval file to be created

-v, --vcf
  Name of the input VCF file
  
-p, --print
  Optional flag to print out list of SVs to stdout
```
## Pileup
```
$ muCNV pileup -s <sample ID> -v <VCF> -V <Interval> -f <GRCh file> -b <BAM/CRAM file> -o <output prefix>

   -o <string>,  --outdir <string>
     Output directory, current directory if omitted

   -s <string>,  --sample <string>
     (required)  Sample ID, also used for output filenames (.pileup, .var,
     .idx)

   -f <string>,  --gcFile <string>
     File containing GC content information (default: GRCh38.gc)

   -V <string>,  --interVal <string>
     Binary interval file containing candidate SVs

   -v <string>,  --vcf <string>
     VCF file containing candidate SVs

   -b <string>,  --bam <string>
     (required)  Input BAM/CRAM file name
```

 - GC content file for human reference genome build 38 is provided in resources/GRCh38.v3.gc
 - Either VCF file (-v) or binary interval file (-V) is required
 
## Merge
```
$ muCNV merge -I [input.list] -o [output_name] -f GRCh38.gc

   -V <string>,  --interVal <string>
     (required)  Binary interval file containing candidate SVs

   -o <string>,  --output <string>
     (required)  Output base filename for merged pileup

   -i <string>,  --index <string>
     (required)  Text file containing list of pileup samples

   -f <string>,  --gcFile <string>
     (required)  File containing GC content information
```

## Genotype
```
$ muCNV genotype [-i <string>] [-f <string>] [-V <string>] [-v <string>] [-o <string>] [-l] [-n <integer-integer>] 

   -p <number(0-1.0)>,  --pmax <number(0-1.0)>
     Maximum overlap between depth clusters

   -x <string>,  --exclude <string>
     List of sample IDs to be excluded from genotyping

   -n <integer-integer>,  --numbers <integer-integer>
     variants in range (from-to) only from the binary SV index file

   -r <chr:startpos-endpos>,  --region <chr:startpos-endpos>
     Genotype specific genomic region

   -c <integer (1-24)>,  --chr <integer (1-24)>
     Chromosome number (1-24) if pileup contains only a single chromosome (default for merged)

   -i <string>,  --index <string>
     Text file containing list of pileups (required)

   -f <string>,  --gcFile <string>
     File containing GC content information

   -V <string>,  --interVal <string>
     Binary interval file containing candidate SVs

   -v <string>,  --vcf <string>
     VCF file containing candidate SVs

   -o <string>,  --out <string>
     Output filename

   -l,  --lessheader
     Do not print header in genoptyed VCF

   -a,  --all
     Include all variants that failed genotyping in the output VCF
```    
 - GC content file for human reference genome build 38 is provided in resources/GRCh38.v3.gc
 - Either VCF file (-v) or binary interval file (-V) is required
 - GRCh38 and Interval file should be the same to the ones used for pileup and merging
