SEQ_FILE=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print}' fq_files/ids` #match SLURM array to list of sequence IDs
REFERENCE=H37Rv.fasta
PICARD=picard.jar

ml bcftools
ml samtools
ml bwa
ml intel/2022.2-dpcpp
ml python/3.12.0

echo $SEQ_FILE

#need to align to reference genome for lineage calling (don't use ancestral reference here)
rgroup="@RG\\tID:"${SEQ_FILE}"\\tSM:"${SEQ_FILE}"_sm\\tPU:"${SEQ_FILE}"_pu\\tLB:"${SEQ_FILE}"_lb"
bwa mem -t 40 -R ${rgroup} ${REFERENCE} ${SEQ_FILE}.P1.filtered.fastq ${SEQ_FILE}.P2.filtered.fastq | \
	awk '$1 ~ /^@/ || $5 == 60' | \
	samtools view -bt ${REFERENCE} - | \
	samtools sort -o ${SEQ_FILE}.H37Rv.bam
	
java -jar ${PICARD} MarkDuplicates I=${SEQ_FILE}.H37Rv.bam O=${SEQ_FILE}.nd.H37Rv.bam M=${SEQ_FILE}.H37Rv.dup.metrix \
	ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

samtools index ${SEQ_FILE}.nd.H37Rv.bam

#fast lineage caller requires a VCF file - need to create one -mv to only include variant positions or -m for all positions
bcftools mpileup -Ov -f ${REFERENCE} ${SEQ_FILE}.nd.H37Rv.bam > ${SEQ_FILE}.fullH37Rv.vcf 

fast-lineage-caller ${SEQ_FILE}.fullH37Rv.vcf --count --out ${SEQ_FILE}.lineage.tsv 


