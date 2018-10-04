#Runs remapping of ENCODE data on CoBrain
 
file_str ='''ENCFF670LIE
ENCFF636QII
ENCFF074BOV
ENCFF486POD
ENCFF887ZOX
ENCFF148RNW
ENCFF890WWJ
ENCFF131JRU
ENCFF280XOG
ENCFF064IOO
ENCFF321CCY
ENCFF126QEZ
'''

from ptes.lib.general import init_file, writeln_to_file, shell_call

bam_folder = '/uge_mnt/home/dp/ngs/encodedcc/hg19/ENCODE/'
segemehl_bin = '/uge_mnt/home/sunnymouse/tools/segemehl_0_2_0/segemehl'        

for i, name in enumerate(file_str.split()):
    id = name.strip(None)                    
    folder_name = '/uge_mnt/home/sunnymouse/projects/PTES/ENCODE/%s' % id   #without end 
    sh_filename = 'segemehl_%i.sh' % i    
    init_file(sh_filename)
    cmd_list = []
    cmd_list.append('#!/bin/bash -il \ncd %s \n' %  '/uge_mnt/home/sunnymouse/projects/PTES/ENCODE/')
    cmd_list.append('source /uge_mnt/home/sunnymouse/tools/miniconda2/bin/activate')
    '''
    cmd_list.append('samtools collate -uOn 256 %s%s.bam %s/tmp-prefix \
                    | samtools fastq - > %s/%s.fq' % (bam_folder, id, folder_name, folder_name, id))
    cmd_list.append('gzip %s/%s.fq' % (folder_name, id))                
    cmd_list.append('%s/segemehl.x -S \
                                -Z 10  \
                                -t 10 \
                                -s \
                                -i /uge_mnt/home/sunnymouse/Human_ref/GRCh37.p13.genome.idx \
                                -d /uge_mnt/home/sunnymouse/Human_ref/GRCh37.p13.genome.fa  \
                                -q %s/%s.fq \
                                -o %s/segemehl.sam \
                                -u %s/segemehl_unmapped' % (segemehl_bin, folder_name,id, folder_name, folder_name)) 
    cmd_list.append('samtools view %s/segemehl.sam > %s/segemehl.sam.nohead' % (folder_name, folder_name))       
    cmd_list.append('samtools view -b %s/segemehl.sam > %s/segemehl.bam' % (folder_name, folder_name))    
    cmd_list.append('%s/segemehl.sam' % folder_name)    
    '''
    cmd_list.append('python segemehl_encode.py -i %s/segemehl.sam.nohead -o %s' % (folder_name, folder_name))    
    
    writeln_to_file('\n'.join(cmd_list), sh_filename)

    shell_call('chmod +x ./%s' % sh_filename)
   