# Detection of WBV reads from WGS
This repo contains a Python script which detects EBV reads from WGS data (BAM) with the chrEBV in the reference human genome.  

## Overview
- This script reads a BAM file mapping to GRCh38 carrying the EBV genome, chrEBV.  
- It extracts reads mapping to chrEBV and outputs below.  
1) number of reads mapping to chrEBV  
2) number of EBV read pairs  
3) number of EBV read pairs mapping outside of BamHI W repeat  
4) number of EBV read pairs mapping inside of BamHI W repeat  
5) BamHI W repeat copy number estimate  
6) total mapped bases of chrEBV.  
- Properly-mapped read pair with mapping length greater than 145 for each read 1 and 2 is considered as a EBV read.
- EBV genome contains a repeat region, called BamHI W repeat. Reads derived from this region can be multimaping due to the repetitive nature. To avoid the counting of multi-mapping, this script extracts reads mapping to BamHI W repeat (positions 12,000–35,355 bp) and counts the number of read pairs accounting for multi-mapping.
- BAM file can contain a very high number of EBV reads due to viremia. In such corner cases, the analysis time will be longer. To avoid the long run time, it downsamples reads if the number of reads mapping to chrEBV is greater than 350,000 reads and analyzes downsampled reads.
- This script also outputs SNVs in chrEBV.


## Usage
```
# Help message
python detect_ebv_reads.py -h

# analyze "sample1.bam" and output results in "sample1.tsv"
python detect_ebv_reads.py \
-i sample1.bam \
-o sample1.tsv \
-c chrEBV
```


## Prerequisites
- Python 3.7 or higher
- samtools (tested with v1.21)
- pysam (tested with v0.22.1)
- pandas (tested with v2.2.2)
- numpy (tested with v1.26.4)


