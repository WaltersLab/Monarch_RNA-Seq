#### RNA-seq reads mapping using STAR ####
# Reference[1] https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
# Reference[2] https://www.youtube.com/watch?v=Ju6PtQD-H34

## Set directories 
dir='/scratch/whtan/RNASeq'
dir_r='/scratch/whtan/Ref'
dir_in='/scratch/whtan/RNASeq/RNASeq_reads/Clean'
dir_out='/scratch/whtan/RNASeq/Outputs'
cd $dir

## Generate genome indexes 
# Monarch genome and annotation files obtained from MonarchBase http://monarchbase.umassmed.edu/home.html
mkdir $dir_r/RNAseq_genome/
STAR --runThreadN 8 --runMode genomeGenerate \
--genomeDir $dir_r/RNAseq_genome/ --genomeFastaFiles $dir_r/Dp_genome_v3.fasta \
--sjdbGTFfile $dir_r/Dp_geneset_OGS2.gtf --sjdbGTFfeatureExon CDS \
--sjdbOverhang 100

## Run mapping for each sample through loop 
# Read a file that contains a list of samlpe names
filename="Sample_list.txt"
exec < $filename
while read sample
do
mkdir $dir_out/$sample
# Align reads to the reference
STAR --runThreadN 8 --runMode alignReads \
--genomeDir $dir_r/RNAseq_genome/ \
--readFilesIn $dir_in/$sample/*_1.fq.gz --readFilesCommand zcat \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix $dir_out/$sample/$sample_
done 