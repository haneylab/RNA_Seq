#!/usr/bin/bash

#load software
module load rsem/1.3.0
module load bowtie2/2.2.9

export rna_home=#Location of RNA-Seq files


#Prepare RSEM reference
rsem-prepare-reference \
--bowtie2 \
$rna_home/ref/ \
$rna_home/ref/A_thal 


#echo Prepare reference done

#TPM caculation

while IFS='' read -r samp || [ -n "$samp" ]; do #loop through the samples.txt file

mkdir $rna_home/${samp}
cd $rna_home/${samp}
rsem-calculate-expression \
-p 8 \
--bowtie2 \
--paired-end \
$rna_home/fastq/CCTC1ANXX_8_1_${samp}_75bp.concat_chastity_passed.fastq \
$rna_home/fastq/CCTC1ANXX_8_2_${samp}_75bp.concat_chastity_passed.fastq \
$rna_home/ref/A_thal \
${samp}

done < $rna_home/samples.txt

