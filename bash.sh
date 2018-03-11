#!/bin/bash	

THREADS=(number_of_threads_available)
DATASET=(name_of_dataset)



# --------------------------------------------------------------------------------------------------------------------------------
# find_circ
# --------------------------------------------------------------------------------------------------------------------------------

bowtie2 -p $THREADS --very-sensitive --mm -M20 --score-min=C,-15,0 -x /path/to/bowtie2_index -q -1 /path/to/file1.fq.gz -2 /path/to/file2.fq.gz | samtools view -hbuS - | samtools sort - output

samtools view -hf 4 output.bam | samtools view -Sb - > unmapped.bam 

python unmapped2anchors.py unmapped.bam | gzip > anchors.qfa.gz

bowtie2 -p $THREADS --reorder --mm -M20 --score-min=C,-15,0 -q -x /path/to/bowtie2_index -U anchors.qfa.gz | python find_circ.py -G /path/to/chomosomes.fa -p prefix -s find_circ.sites.log > find_circ.sites.bed 2> find_circ.sites.reads

grep circ find_circ.sites.bed | grep -v chrM | python sum.py -2,3 | python scorethresh.py -16 1 | python scorethresh.py -15 2 | python scorethresh.py -14 2 | python scorethresh.py 7 2 | python scorethresh.py 8,9 35 | python scorethresh.py -17 100000 > $DATASET.findcirc.bed


#for find_circ(40x40):

grep circ find_circ.sites.bed | grep -v chrM | python sum.py -2,3 | python scorethresh.py -16 1 | python scorethresh.py -15 2 | python scorethresh.py -14 2 | python scorethresh.py 7 2 | python scorethresh.py 8 40 | python scorethresh.py 9 40 | python scorethresh.py -17 100000 > $DATASET.findcirc_40x40.bed



# --------------------------------------------------------------------------------------------------------------------------------
# circRNA_finder
# --------------------------------------------------------------------------------------------------------------------------------

STAR --genomeDir /path/to/star_index \
     --readFilesIn /path/to/file1.fq /path/to/file2.fq \
     --runThreadN $THREADS \
     --chimSegmentMin 20 \
     --chimScoreMin 1 \
     --alignIntronMax 100000 \
     --outFilterMismatchNmax 4 \
     --alignTranscriptsPerReadNmax 100000 \
     --outFilterMultimapNmax 2 \
     --outFileNamePrefix /path/to/star/output


perl postProcessStarAlignment.pl /path/to/star/output /path/to/circRNA_finder/output



# --------------------------------------------------------------------------------------------------------------------------------
# CIRI
# --------------------------------------------------------------------------------------------------------------------------------

bwa mem -t $THREADS -T 19 -v 2 /path/to/bwa_index /path/to/file1.fq /path/to/file2.fq >output.sam

perl CIRI_v1.2.pl -I output.sam -O /path/to/ciri/output.txt -F /path/to/genome.fa -P



# --------------------------------------------------------------------------------------------------------------------------------
# CIRI2
# --------------------------------------------------------------------------------------------------------------------------------

bwa mem -t $THREADS -T 19 -v 2 /path/to/bwa_index /path/to/file1.fq /path/to/file2.fq >output.sam

perl CIRI_v2.0.6.pl -T $THREADS -I output.sam -O /path/to/ciri/output.txt -F /path/to/genome.fa


# --------------------------------------------------------------------------------------------------------------------------------
# circExplorer
# --------------------------------------------------------------------------------------------------------------------------------

tophat2 -a 6 --microexon-search -m 2 -p $THREADS -o circExplorer_output /path/to/bowtie1_index /path/to/file1.fq.gz /path/to/file2.fq.gz

bamToFastq -i circExplorer_output/unmapped.bam -fq circExplorer_output/unmapped.fastq

tophat2 -o circExplorer_output -p $THREADS --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search path/to/bowtie2_index circExplorer_output/unmapped.fastq

python CIRCexplorer.py -f circExplorer_output/accepted_hits.bam -g /path/to/genome.fa -r /path/to/gene_annotation_ucsc_hg19.txt -o output.txt 


# --------------------------------------------------------------------------------------------------------------------------------
# circExplorer2
# --------------------------------------------------------------------------------------------------------------------------------


python command_parse.py align -p $THREADS -o circExplorer2_output -G /path/to/gtf -i /path/to/genome_index -j /path/to/genome_index /path/to/file1.fq.gz /path/to/file2.fq.gz
python command_parse.py parse -t TopHat-Fusion /path/to/tophat_fusion/accepted_hits.bam
python command_parse.py annotate -r /path/to/gene_annotation_ucsc_hg19.txt -g /path/to/genome.fa /path/to/output_folder
python command_parse.py assemble -p $THREADS -r /path/to/gene_annotation_ucsc_hg19.txt /path/to/output_folder
python command_parse.py denovo -r /path/to/gene_annotation_ucsc_hg19.txt -g /path/to/genome.fa /path/to/output_folder

cat /path/to/output_folder/denovo/*_circ.txt >> /path/to/output_folder/$DATASET.CIRCexplorer2.circ.txt


# --------------------------------------------------------------------------------------------------------------------------------
# MapSplice
# --------------------------------------------------------------------------------------------------------------------------------

python MapSplice-v2.1.8/mapsplice.py -1 /path/to/file1.fq -2 file2.fq -c /path/to/chromosome_sequences -x /path/to/bowtie1_index -p $THREADS -o Mapsplice_output --min-fusion-distance 200 --gene-gtf Homo_sapiens.GRCh37.66.gtf --fusion 



# --------------------------------------------------------------------------------------------------------------------------------
# KNIFE
# --------------------------------------------------------------------------------------------------------------------------------


sh findCircularRNA_SLURM.sh /path/to/srr complete /path/to/output $DATASET 8 sam_large_phred33
sh findCircularRNA_SLURM.sh /path/to/srr complete /path/to/output $DATASET 8 sam_large_phred33_unaligned
sh findCircularRNA_SLURM.sh /path/to/srr complete /path/to/output $DATASET 8 sam_large_phred33_R2analysis
python analysis/combineSwappedReadsGLM.py -a /path/to/output/$DATASET/circReads -b /path/to/output/$DATASETSwapped/circReads -q complete -v



# creating mock annotation reference

python makeExonDB.py -f /path/to/genome.fasta -a /path/to/mock.gtf -o /path/to/db_output
sh createJunctionIndex.sh /path/to/KNIFE /path/to/output mock 1000000


# --------------------------------------------------------------------------------------------------------------------------------
# DCC
# --------------------------------------------------------------------------------------------------------------------------------

mkdir sample1
cd sample1
STAR --runThreadN $THREADS \
       --genomeDir /path/to/star_index \
       --outSAMtype BAM SortedByCoordinate \
       --readFilesIn /path/to/file1.fq.gz /path/to/file2.fq.gz \
       --readFilesCommand zcat \
       --outFileNamePrefix dataset_name \
       --outReadsUnmapped Fastx \
       --outSJfilterOverhangMin 15 15 15 15 \
       --alignSJoverhangMin 15 \
       --alignSJDBoverhangMin 15 \
       --outFilterMultimapNmax 20 \
       --outFilterScoreMin 1 \
       --outFilterMatchNmin 1 \
       --outFilterMismatchNmax 2 \
       --chimSegmentMin 15 \
       --chimScoreMin 15 \
       --chimScoreSeparation 10 \
       --chimJunctionOverhangMin 15 \

cd ..
mkdir mate1
cd mate1
STAR --runThreadN $THREADS \
       --genomeDir /path/to/star_index \
       --outSAMtype None \
       --readFilesIn /path/to/file1.fq.gz \
       --readFilesCommand zcat \
       --outFileNamePrefix dataset_name \
       --outReadsUnmapped Fastx \
       --outSJfilterOverhangMin 15 15 15 15 \
       --alignSJoverhangMin 15 \
       --alignSJDBoverhangMin 15 \
       --seedSearchStartLmax 30 \
       --outFilterMultimapNmax 20 \
       --outFilterScoreMin 1 \
       --outFilterMatchNmin 1 \
       --outFilterMismatchNmax 2 \
       --chimSegmentMin 15 \
       --chimScoreMin 15 \
        --chimScoreSeparation 10 \
       --chimJunctionOverhangMin 15 \

cd ..

mkdir mate2
cd mate2
STAR --runThreadN $THREADS \
       --genomeDir $GENOME \
       --outSAMtype None \
       --readFilesIn $PATH2 \
       --readFilesCommand zcat \
       --outFileNamePrefix $DATASET \
       --outReadsUnmapped Fastx \
       --outSJfilterOverhangMin 15 15 15 15 \
       --alignSJoverhangMin 15 \
       --alignSJDBoverhangMin 15 \
       --seedSearchStartLmax 30 \
       --outFilterMultimapNmax 20 \
       --outFilterScoreMin 1 \
       --outFilterMatchNmin 1 \
       --outFilterMismatchNmax 2 \
       --chimSegmentMin 15 \
       --chimScoreMin 15 \
       --chimScoreSeparation 10 \
       --chimJunctionOverhangMin 15 \

cd ..

COJ="Chimeric.out.junction"

echo "/path/to/sample1/$DATASET$COJ" > samplesheet
echo "/path/to/mate1/$DATASET$COJ" > mate1file
echo "/path/to/mate2/$DATASET$COJ" > mate2file

python /home/tbhmbi/DCC/DCC/main.py @samplesheet -mt1 @mate1file -mt2 @mate2file -D -Pi -F -M -Nr 1 1 -fg -A $FASTA -N -T $THREADS



# --------------------------------------------------------------------------------------------------------------------------------
# ACSF
# --------------------------------------------------------------------------------------------------------------------------------


perl change_fastq_header.pl /path/to/file1.fq file1_changed.fa Truseq_file1left
perl change_fastq_header.pl /path/to/file2.fq  file2_changed.fa Truseq_file2left
perl Truseq_merge_unique_fa.pl UNMAP dataset_name file1_changed.fa file2_changed.fa

echo "BWA_folder	/path/to/bwa" > SPEC.txt
echo "BWA_genome_Index	/path/to/bwa/index/genome.fa" >> SPEC.txt
echo "BWA_genome_folder	/path/to/bwa/index" >> SPEC.txt
echo "ACF_folder	/path/to/acfs" >> SPEC.txt
echo "CBR_folder	/path/to/CB_splice/" >> SPEC.txt
echo "Agtf	/path/to/gtf/" >> SPEC.txt
echo "UNMAP	UNMAP" >> SPEC.txt
echo "UNMAP_expr	UNMAP_expr" >> SPEC.txt
echo "Seq_len	100" >> SPEC.txt
echo "Thread	$THREADS" >> SPEC.txt
echo "Strandness	no" >> SPEC.txt

perl ACF_MAKE.pl SPEC.txt BASH_file.sh

chmod +x BASH_file.sh

bash BASH_file.sh 



# --------------------------------------------------------------------------------------------------------------------------------
# UROBORUS
# --------------------------------------------------------------------------------------------------------------------------------


tophat -p $THREADS -o /path/to/output /path/to/bowtie2_index $PATH1 $PATH2

samtools view /path/to/output/unmapped.bam > /path/to/output/unmapped.sam

perl UROBORUS_0.1.2.pl -temp -p $THREADS -index /path/to/bowtie1_index -gtf /path/to/gtf -fasta /path/to/fastafiles /path/to/output/unmapped.sam


