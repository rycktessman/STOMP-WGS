ml bcftools/1.15.1
ml vcftools

NUM_ISOLATES=241

#merge vcfs
bcftools merge *vcf.gz --force-samples -0 -O v -o merged.vcf

#filter and remove indels:
vcftools --vcf merged.vcf --exclude-bed ${BIN}/filter_tree_bed.txt --remove-indels --recode --recode-INFO-all --out remove_indel

#transform to phylip:
ml intel/2022.2-dpcpp
ml python/3.12.0
python vcf2phylip.py --input remove_indel.recode.vcf --output-prefix remove_indel -m ${NUM_ISOLATES}

# tree building
iqtree-1.6.12-Linux/bin/iqtree -s remove_indel.min${NUM_ISOLATES}.phy -m MFP -mset GTR -alrt 1000 
