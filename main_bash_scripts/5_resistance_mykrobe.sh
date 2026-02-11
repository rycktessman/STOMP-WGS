module load python/3.8.6
source PATH/myenv/bin/activate #if needed, activate python environment that contains mykrobe installation

SEQ_FILE=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print}' fq_files/ids`  #match SLURM array to list of sequence IDs

mykrobe predict --sample $SEQ_FILE --species "tb" -1 $SEQ_FILE.nd.H37Rv.bam  --output $SEQ_FILE.mykrobe.csv 
