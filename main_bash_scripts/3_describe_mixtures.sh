module purge
ml kraken2
ml ncurses
ml gcc/9.3.0
ml intel/2022.2-dpcpp
ml python/3.12.0

VCFMIX=pipeline_files/VCFMIX #location of VCFMIX

SEQ_FILE=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print}' fq_files/ids` #match SLURM array to list of sequence IDs
echo ${SEQ_FILE}

python describe_mixtures.py -p ${VCFMIX} -i ${SEQ_FILE} -f ${SEQ_FILE}

