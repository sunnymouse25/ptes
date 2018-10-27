#Runs remapping of ENCODE data on CoBrain
 
file_str ='''ENCFF670LIE
ENCFF636QII
ENCFF074BOV
ENCFF486POD
ENCFF887ZOX
ENCFF148RNW
ENCFF890WWJ
ENCFF280XOG
ENCFF064IOO
ENCFF321CCY
ENCFF126QEZ
'''

from ptes.lib.general import init_file, writeln_to_file, shell_call

for i, name in enumerate(file_str.split()):
    id = name.strip(None)                    
    folder_name = '/home/sunnymouse/projects/PTES/ENCODE/%s' % id   #without end 
    sh_filename = 'star_%i.sh' % i    
    init_file(sh_filename)
    cmd_list = []
    '''
    cmd_list.append('#!/bin/bash -il \ncd %s \n' %  '/home/sunnymouse/projects/PTES/ENCODE/')
    cmd_list.append('source /home/sunnymouse/miniconda22/bin/activate')    
    '''
    cmd_list.append('STAR   \
    --runThreadN 10   \
    --genomeDir /home/sunnymouse/Human_ref/hg19_chrom/Genome_indices/   \
    --genomeLoad NoSharedMemory   \
    --readFilesIn %s/%s.fq1.gz\
    --readFilesCommand gunzip -c \
    --outFileNamePrefix %s/mate1_  \
    --outSAMtype SAM      \
    --outSAMattributes NH   HI   AS   NM   MD      \
    --outSAMunmapped None      \
    --outSAMheaderHD @HD   VN:1.4   SO:coordinate      \
    --outFilterType BySJout   \
    --outFilterMultimapNmax 20   \
    --outFilterMismatchNmax 999   \
    --outFilterMismatchNoverReadLmax 0.04   \
    --alignIntronMax 1000000   \
    --alignMatesGapMax 1000000   \
    --alignSJoverhangMin 8   \
    --alignSJDBoverhangMin 1   \
    --chimSegmentMin 10   \
    --sjdbGTFfile /home/sunnymouse/Human_ref/gencode.v19.annotation.gtf   \
    --sjdbScore 3 \
    --chimMainSegmentMultNmax 1 \
    --chimJunctionOverhangMin 5 \
    --chimMultimapNmax 0 \
    --twopassMode Basic' % (folder_name, id, folder_name))
    cmd_list.append('STAR   \
    --runThreadN 10   \
    --genomeDir /home/sunnymouse/Human_ref/hg19_chrom/Genome_indices/   \
    --genomeLoad NoSharedMemory   \
    --readFilesIn %s/%s.fq2.gz\
    --readFilesCommand gunzip -c \
    --outFileNamePrefix %s/mate2_  \
    --outSAMtype SAM      \
    --outSAMattributes NH   HI   AS   NM   MD      \
    --outSAMunmapped None      \
    --outSAMheaderHD @HD   VN:1.4   SO:coordinate      \
    --outFilterType BySJout   \
    --outFilterMultimapNmax 20   \
    --outFilterMismatchNmax 999   \
    --outFilterMismatchNoverReadLmax 0.04   \
    --alignIntronMax 1000000   \
    --alignMatesGapMax 1000000   \
    --alignSJoverhangMin 8   \
    --alignSJDBoverhangMin 1   \
    --chimSegmentMin 10   \
    --sjdbGTFfile /home/sunnymouse/Human_ref/gencode.v19.annotation.gtf   \
    --sjdbScore 3 \
    --chimMainSegmentMultNmax 1 \
    --chimJunctionOverhangMin 5 \
    --chimMultimapNmax 0 \
    --twopassMode Basic' % (folder_name, id, folder_name))
    cmd_list.append("cat %s/mate1_Chimeric.out.junction | awk '$1 ==$4 && $3 ==$6 && $7 >= 0 && $8+$9<=5 ' | sort |  uniq  > %s/mate1_Chimeric.out.junction.filtered" % (folder_name, folder_name))                
    cmd_list.append("cat %s/mate2_Chimeric.out.junction | awk '$1 ==$4 && $3 ==$6 && $7 >= 0 && $8+$9<=5 ' | sort |  uniq  > %s/mate2_Chimeric.out.junction.filtered" % (folder_name, folder_name))                   
    cmd_list.append('mkdir %s/mate1' % folder_name)    
    cmd_list.append('mkdir %s/mate2' % folder_name)    
    cmd_list.append('python star_encode_SE.py \
                    -i %s/mate1_Chimeric.out.junction.filtered\
                    -s %s/mate1_Aligned.out.sam \
                    -o %s/mate1\
                    -g /home/sunnymouse/Human_ref/GRCh37.p13.genome.fa\
                    -gtf /home/sunnymouse/Human_ref/hg19_exons.gtf' % (folder_name, folder_name, folder_name))    
    cmd_list.append('python star_encode_SE.py \
                    -i %s/mate2_Chimeric.out.junction.filtered \
                    -s %s/mate2_Aligned.out.sam.sam \
                    -o %s/mate2\
                    -g /home/sunnymouse/Human_ref/GRCh37.p13.genome.fa\
                    -gtf /home/sunnymouse/Human_ref/hg19_exons.gtf' % (folder_name, folder_name, folder_name))    
        
    
    writeln_to_file('\n'.join(cmd_list), sh_filename)

    shell_call('chmod +x ./%s' % sh_filename)
   