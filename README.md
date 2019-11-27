# muCNV

Multi-sample SV genotyper for large-scale WGS data

Version 0.9.5
Last edited: Nov. 27, 2019, (c) Goo Jun

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
pileup
merge
genotype
print
vcf2int
gcidx
filter
```    

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
$ muCNV pileup -s <sample ID> -v <VCF> -f <GRCh file> -b <BAM/CRAM file> -o <output prefix>

-o, --outdir
  Output directory. Default is current (.) directory
  
-s,  --sample
  Sample ID for output filename base

-f <string>,  --gcFile <string> 
  File containing GC content information

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

-o <string>,  --output <string> 
  (required)  Output base filename for merged pileup

-i <string>,  --index <string>
  (required)  Text file containing list of pileup samples

-f <string>,  --gcFile <string>
  File containing GC content information

-V <string>,  --interVal <string>
  Binary interval file containing candidate SVs
```

## Genotype
```
$ muCNV genotype [-i <string>] [-f <string>] [-V <string>] [-v <string>] [-o <string>] [-l] [-n <integer-integer>] 

-r <chr:startpos-endpos>,  --region <chr:startpos-endpos>
  Genotype specific genomic region

-n <integer-integer>,  --numbers <integer-integer>
  variants in range (from-to) only

-i <string>,  --index <string>
   List file containing list of intermediate pileups. Required.

-f <string>,  --gcFile <string>
  File containing GC content information. Required.

-V <string>,  --interVal <string>
  Binary interval file containing candidate SVs
   
-v <string>,  --vcf <string>
  VCF file containing candidate SVs

-p <number>, --pmax <number>
  Maximum overlap between depth cluster PDFs

-r <chr:startpos-endpos>,  --region <chr:startpos-endpos>
  Genotype specific genomic region

-o <string>,  --out <string>
  Output filename

-l,  --lessheader
  Do not print header in genoptyed VCF
```    
 - GC content file for human reference genome build 38 is provided in resources/GRCh38.v3.gc
 - Either VCF file (-v) or binary interval file (-V) is required
 - GRCh38 and Interval file should be the same to the ones used for pileup and merging
