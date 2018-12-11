# muCNV

Multi-sample SV genotyper for large-scale WGS data

Last edited: Dec. 11, 2018, (c) Goo Jun

muCNV has two-stages for multi-sample SV genotyping:
1. Generate *pileups* from CRAM/BAM files using candidate list of SVs
2. Multi-sample genotyping from summary files.

There is also an optional 'pileup merging' step in-between, which is recommended for >300 samples.

## Basis usage:

    $ mucnv [command] [options]
    $ mucnv [command] --help
    
 
  

List of commands:

    pileup
    merge
    genotype
    print
    vcf2int
    gcidx
    filter
    

## Pileup

     $ mucnv pileup -s <sample ID> -v <VCF> -f <GRCh file> -b <BAM/CRAM file> -o <output prefix>
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
     

 - GC content file for human reference genome build 38 is provided in resources/GRCh38.gc
 - Either VCF file (-v) or binary interval file (-V) is required
 
 
  
## Merge
  

    $ muCNV merge -I [input.list] -o [output_name] -f GRCh38.gc
    -o <string>,  --output <string> (required)  Output base filename for merged pileup
    -i <string>,  --index <string>
    (required)  Text file containing list of pileup samples
    -f <string>,  --gcFile <string>
    File containing GC content information

## Genotype

    $ muCNV genotype [-i <string>] [-f <string>] [-V <string>] [-v <string>] [-o <string>] [-l] [-n <integer-integer>] 
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
    -o <string>,  --out <string>
    Output filename
    -l,  --lessheader
    Do not print header in genoptyed VCF
    

 - GC content file for human reference genome build 38 is provided in resources/GRCh38.gc
 - Either VCF file (-v) or binary interval file (-V) is required
 - GRCh38 and Interval file should be the same to the ones used for pileup and merging
