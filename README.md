# muCNV

Multi-sample SV genotyper for large-scale WGS data

Version 0.9.7 (TOPMed Freeze 8 SV version)
Last edited: June 28, 2019, (c) Goo Jun

muCNV uses multiple steps for multi-sample SV genotyping, to handle large number of samples and to enable efficient parallelization:
0. Generate binary interval file from input VCF 
1. Generate *pileups* from CRAM/BAM files using candidate list of SVs
2. Merge individual pileups to generate 'merged' pileups
3. Multi-sample genotyping from merged pileups

* The 'merging' step is optional, but is recommended for ~ >500 samples. 
* By default, merging will create pileups separated by chromosomes. 
* When genotyping merged pileups, do not add 'chr' postfixes in your 'list' file

## Basis usage:
```
$ mucnv [command] [options]
$ mucnv [command] --help
```    
List of commands:
```
vcf2int
pileup
merge
genotype
gcidx
```    

##   Generating variant list from input VCF
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
$ muCNV pileup -s <sample ID> -v <VCF> -f <GRCh file> -b <BAM/CRAM file> -o <output prefix>

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
