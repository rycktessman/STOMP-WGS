module purge
ml kraken2
ml ncurses
ml gcc/9.3.0
ml intel/2022.2-dpcpp
ml python/3.12.0

REFERENCE="MTB_ancestor_reference.fna"
KRAKENDB_KEY="27042024"
KRAKENDB="kraken-db-27042024"
SEQTK="seqtk-1.2/seqtk"
BWA="bwa_0713/bwa"
PICARD=picard.jar
VARSCAN=VarScan.v2.3.7.jar
BEDTOOLS="bedtools2-2.26.0/bin/genomeCoverageBed"
SAMTOOLS="samtools" #version 1.3.1
GATK="gatk-4.0.2.1/gatk-package-4.0.2.1-local.jar"
ANNOTATION="annotation/H37Rv.annotation.tab"

TYPE=1 #1 or 2, for 2 possible filename/ID formats \\

if["$TYPE"==1]; then
	IDs=$(ls *_R1.fastq.gz | egrep '_R1.fastq.gz' | sed s/'_R1.fastq.gz'//)
fi

if["$TYPE"==2]; then
	IDs=$(ls *_R1_001.fastq.gz | egrep '_R1_001.fastq.gz' | sed s/'_R1_001.fastq.gz'//)
fi


#loop over samples
for SAMPLE in IDs
do
	echo ${SAMPLE}

	#trim the reads 
	if [ "$TYPE"==2]; then
		java -jar trimmomatic-0.36.jar PE -threads 30 -phred33 ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz ${SAMPLE}.P1.clean.fastq \
			${SAMPLE}.U1.clean.fastq ${SAMPLE}.P2.clean.fastq ${SAMPLE}.U2.clean.fastq \
			TRAILING:10 SLIDINGWINDOW:25:20 MINLEN:50 2> ${SAMPLE}.trimmomatic.log
	fi

	if [ "$TYPE"==1 ]; then
		java -jar trimmomatic-0.36.jar PE -threads 30 -phred33 ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz ${SAMPLE}.P1.clean.fastq \
			${SAMPLE}.U1.clean.fastq ${SAMPLE}.P2.clean.fastq ${SAMPLE}.U2.clean.fastq \
			TRAILING:10 SLIDINGWINDOW:25:20 MINLEN:50 2> ${SAMPLE}.trimmomatic.log
	fi

	#map reads to kraken database
	kraken2 --db ${KRAKENDB} --use-names --report ${SAMPLE}.${KRAKENDB_KEY}.kraken.report \
		--paired ${SAMPLE}.P1.clean.fastq ${SAMPLE}.P2.clean.fastq > ${SAMPLE}.kraken.${KRAKENDB_KEY}

	#extract names of Mtb reads from kraken output
	python  extract_kraken_reads.py \
		-k ${SAMPLE}.kraken.${KRAKENDB_KEY} -s ${SAMPLE}.P1.clean.fastq -s2 ${SAMPLE}.P2.clean.fastq \
		-r ${SAMPLE}.${KRAKENDB_KEY}.kraken.report -t 77643 --include-children --fastq-output \
		-o ${SAMPLE}.P1.filtered.fastq -o2 ${SAMPLE}.P2.filtered.fastq > ${SAMPLE}.krakenpy.log

	#map to reference genome, generate SAM file
	rgroup="@RG\\tID:"${SAMPLE}"\\tSM:"${SAMPLE}"_sm\\tPU:"${SAMPLE}"_pu\\tLB:"${SAMPLE}"_lb"
	${BWA} mem -t 40 -R ${rgroup} ${REFERENCE} ${SAMPLE}.P1.filtered.fastq ${SAMPLE}.P2.filtered.fastq | \
		awk '$1 ~ /^@/ || $5 == 60' | \
		${SAMTOOLS}	view -bt ${REFERENCE} - | \
		${SAMTOOLS} sort -o ${SAMPLE}.sort.bam

	#use MarkDuplicates (a tool in Picard) to remove duplicates resulting from non-random fragmentation of the genome
	java -jar ${PICARD} MarkDuplicates I=${SAMPLE}.sort.bam O=${SAMPLE}.nd.sort.bam M=${SAMPLE}.dup.metrix \
		ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

	#Calculate coverage file (# reads by position)
	${BEDTOOLS} -ibam ${SAMPLE}.nd.sort.bam -d -g ${REFERENCE} > ${SAMPLE}.coverage

	#generate meancov (coverage) file and move files out of main folder if they don't pass QC
	python Coverage.py --p ${SAMPLE} --r ${REFERENCE} --f 1

	#only generate mpileup file and SNP calls if sequence passes QC (based on whether the file has been moved or not in previous step)
	#note: QC for % Mtb (from kraken) and for mixed infections is done later (manually exclude those isolates later).
	if [ -e ${SAMPLE}.nd.sort.bam ]; then
		echo "File passed QC"
		#SAMTOOLS generates mpileup file from BAM file
		${SAMTOOLS} index ${SAMPLE}.nd.sort.bam
		${SAMTOOLS} mpileup -AB -f ${REFERENCE} ${SAMPLE}.nd.sort.bam > ${SAMPLE}.mpileup

		#raw variant calling
		java -jar ${VARSCAN} pileup2snp ${SAMPLE}.mpileup --p-value 0.01 --min-coverage 3 --min-reads2 3 --min-freq-for-hom 0.9 --min-var-freq 0.05 > ${SAMPLE}.snp

		#generate vSNPs - they will be used for rescue SNPs
		java -jar ${VARSCAN} pileup2snp ${SAMPLE}.mpileup --p-value 0.01 --min-coverage 10 --min-reads2 6 --min-avg-qual 25 --min-strands2 2 --min-var-freq 0.1 > ${SAMPLE}.vSNPs

		#filtered variants
		java -jar ${VARSCAN} pileup2snp ${SAMPLE}.mpileup --p-value 0.01 --min-coverage 20 --min-reads2 20 --min-avg-qual 25 --min-strands2 2 --min-var-freq 0.9 > ${SAMPLE}.fSNPs 
		
		#indel calling
		java -jar ${VARSCAN} pileup2indel ${SAMPLE}.mpileup --p-value 0.01 --min-freq-for-hom 0.9 --min-var-freq 0.1 > ${SAMPLE}.vscan.indel				
		
		#Variant filtration - annotation filter - run first on vSNPs, then fSNPs, then tSNPs. Returns .annoF files
		python Anno_Filter.py --anno ${ANNOTATION} --snp ${SAMPLE}.vSNPs
		python Anno_Filter.py --anno ${ANNOTATION} --snp ${SAMPLE}.fSNPs
		
		#Variant filtration - indel filter - run first on vSNPs, then fSNPs, then tSNPs. Returns .annoF.delF files
		python Del_Filter.py --snp ${SAMPLE}.vSNPs.annoF --delfile ${SAMPLE}.vscan.indel
		python Del_Filter.py --snp ${SAMPLE}.fSNPs.annoF --delfile ${SAMPLE}.vscan.indel

		#Variant filtration - density filter (window-based filtering) - run only on fSNPs and tSNPs
		python /scratch4/ddowdy1/tess/stomp/code/Density_Filter.py --snp ${SAMPLE}.fSNPs.annoF.delF

		
	else
		echo "File failed QC"
	fi

	touch ${SAMPLE}".done"
done