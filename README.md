# PTES - Post-Transcriptional Exon Shuffling - search tools

## How to install

Use pip to install:
```
pip install git+https://github.com/sunnymouse25/ptes@master --upgrade
 
```

```
pip install git+https://github.com/sunnymouse25/ptes@dev --upgrade
```

## How to use
### Step 1. Obtaining FASTQ files 

You can start either from FASTQ or BAM files. If you have BAM files, you have to convert them to FASTQ with the following command (it is the fastest one):
```
samtools collate -uOn 256 {input.bam} tmp-prefix | samtools fastq -1 {output.fq1} -2 {output.fq2} -
```
### Step 2. Mapping FASTQ files

To obtain chimeric reads in both "mate_inside" and "mate_outside" configurations, you should map your FASTQ files separately, in single-end mode. Use prefixes mate1_ and mate2_ for two FASTQ files, respectively. 
Here is my set of parameters for STAR to map reads, based on ENCODE parameters. I used STAR 2.6.0. Do not forget to obtain STAR genome indices if you don't have them!

```
STAR  
--runThreadN 10  
--genomeDir {STAR_Genome_indices}   
--genomeLoad NoSharedMemory   
--readFilesIn {id.fq1.gz}    
--readFilesCommand gunzip -c 
--outFileNamePrefix mate1_{id}_  
--outSAMtype SAM      
--outSAMattributes NH   HI   AS   NM   MD      
--outSAMunmapped None      
--outSAMheaderHD @HD   VN:1.4   SO:coordinate      
--outFilterType BySJout   
--outFilterMultimapNmax 20   
--outFilterMismatchNmax 999   
--outFilterMismatchNoverReadLmax 0.04   
--alignIntronMax 1000000   
--alignMatesGapMax 1000000   
--alignSJoverhangMin 8  
--alignSJDBoverhangMin 1  
--chimSegmentMin 10   
--sjdbGTFfile {hg19.gtf}   
--sjdbScore 3 
--chimMainSegmentMultNmax 1 
--chimJunctionOverhangMin 5 
--chimMultimapNmax 0 
--chimOutType Junctions 
--twopassMode Basic
```   

### Step 3. Creating lists of files

At this step, for each sample you have 4 files:

```
mate1_Aligned.out.sam
mate2_Aligned.out.sam
mate1_Chimeric.out.junction
mate2_Chimeric.out.junction
```

You have to make 3 files: 2 lists with full paths to these files, one with chimeric files and the other one with SAM files. 
**Important:** similar lines of these two files should contain paths for the **same** sample, but **different** mates. Example:

```
# list_chimeric.txt
/path_to_file/mate1_ENCFF064IOO.Chimeric.out.junction
/path_to_file/mate2_ENCFF064IOO.Chimeric.out.junction
```  

```
# list_sam.txt
/path_to_file/mate2_ENCFF064IOO.Aligned.out.sam
/path_to_file/mate1_ENCFF064IOO.Aligned.out.sam
```  

The third file is list of tags that will be used for grouping results. 
**Important:** similar lines of this file and two previous files should contain paths for the **same** sample. If you have paired-end reads and two FASTQ files, repeat line with tag twice. Example:
```
# list_tag.txt
ENCFF064IOO
ENCFF064IOO
```  

### Step 4. Creating tables with "mate_inside" and "mate_outside"

Run the script **ptes.star_SE_chimeric** to obtain summary tables: the first for each read, the second for each junction. 
   
 
Example:  

```
ptes.star_SE_chimeric
-i list_chimeric.txt
-s list_sam.txt
-t list_tag.txt
-o output/   # output folder
-l 1   # tell the script that the input is lists
-gz 1  # tell the script to gzip the output, optional
-gtf {hg19.gtf}   # full path to genome annotation 
```

Also you can run: 
```
ptes.star_SE_chimeric -h
```
to get the full list of parameters. 

Table for each read (**chim_reads.csv**): each row has the information about one pair of reads mapped in "mate_inside" or "mate_outside" configuration. Table for each junction (**chim_junctions.csv**) groups the first table by coordinates of chimeric junction, each row has the counts of reads mapped in "mate_inside" and "mate_outside" configuration. A JSON file (**junc_dict.json**) has the information about genomic intervals of reads to create BED files.

### Step 5. Creating BED files

Run the script **ptes.junctions_to_bed** to get BED and, optionally, bigBED files for these reads or subset of reads. 

Example to get only reads with annotated junctions:
```bash
ptes.junctions_to_bed 
-t chim_reads.csv   # don't change
-j junc_dict.json.gz # change only to "junc_dict.json" if you disabled gzip option in the previous step
-gz 1 # remove only if you disabled gzip option in the previous step
-q 'annot_donor+annot_acceptor >=1 & letters_ss == "GT/AG" & chim_dist < 30000 & mate_dist < 30000 & mate_dist > -30000' # pandas query to get the subset from chim_reads.csv
-o annot # name of the output folder
-p annot # name of the output prefix for BED files
-sort 1   # sort the output BED files, optional
-bb 1  # make bigBED from the output BED files, optional. bedToBigBed from UCSC should be in your PATH
```

Also you can run: 
```
ptes.junctions_to_bed -h
```
to get the full list of parameters. 

You will get 4 BED files with different number of parts per row:

```bash
$ wc -l annot/
  31928 annot/annot.bed  # all reads, one pair of mates - two rows
  10642 annot/annot.single.bed  # one random read per junction, one pair of mates - two rows
   5922 annot/annot.single.unique.bed  # one random read per junction, one pair of mates - one row
  17769 annot/annot.unique.bed  # all reads, one pair of mates - one row

```    