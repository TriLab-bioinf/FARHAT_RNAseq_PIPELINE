#!/usr/bin/bash
#SBATCH --job-name=RNAseq  # Job name
#SBATCH --output=RNAseq_output_%j.log  # Standard output and error log
#SBATCH --ntasks=1  # Run a single task
#SBATCH --cpus-per-task=16  # Number of CPU cores per task
#SBATCH --mem=64G  # Memory allocation
#SBATCH --time=48:00:00  # Time limit (adjust as needed)
#SBATCH --mail-user=farhat.parween@nih.gov  # Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # Notify on job start, end, or failure

set -o errexit

# Load necessary modules
module load fastqc multiqc trimmomatic STAR samtools picard subread

Help()
{
   echo "This script runs STAR on a fasta file to generate STAR index files."
   echo "Syntax: ${0} [-p|-c|-h]"
   echo "Options:"
   echo "-p     seq file prefix. [./samples.prefix]"
   echo "-h     Prints this help."
}

# Parse command-line arguments
while getopts "hp:d:t:R:a:" option; do
   case $option in
        h) Help; exit;;
        p) input=${OPTARG};;
        c) CONFIG_FILE=${OPTARG};;
        \?) echo "Invalid option"; Help; exit;;
   esac
done

# Set defaults if not provided
input=${input:-./samples.prefix}
CONFIG_FILE=${CONFIG_FILE:-./config.txt}


# Set working directory
WD=$(pwd -P)

# Load configuration variables values from config.txt file
source $CONFIG_FILE

# Run FASTQC
while IFS= read -r prefix; do
    file1=${READS}/${prefix}.R1.fastq.gz
    file2=${READS}/${prefix}.R2.fastq.gz
    fastqc -t ${CPU} $file1 $file2
done < "$input"

multiqc ${READS}/

# Set trimmomatic output directory
TRIM_DIR="${WD}/1-trimmed"
mkdir -p "$TRIM_DIR"

# Trimming with Trimmomatic
while IFS= read -r prefix; do
    file1=${READS}/${prefix}.R1.fastq.gz
    file2=${READS}/${prefix}.R2.fastq.gz
    outP1=${TRIM_DIR}/${prefix}.R1.paired.fastq.gz
    outP2=${TRIM_DIR}/${prefix}.R2.paired.fastq.gz
    outUP1=${TRIM_DIR}/${prefix}.R1.unpaired.fastq.gz
    outUP2=${TRIM_DIR}/${prefix}.R2.unpaired.fastq.gz
    log=${TRIM_DIR}/${prefix}.trim.log
    java -jar  ${TRIMMOMATIC_JAR} PE -threads ${CPU} -trimlog $log $file1 $file2 \
        $outP1 $outUP1 $outP2 $outUP2 ILLUMINACLIP:${TRIMMOMATIC_JARPATH}/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
done < "$input"

# Run FASTQC again on trimmed reads
while IFS= read -r prefix; do
    file1=${TRIM_DIR}/${prefix}.R1.paired.fastq.gz
    file2=${TRIM_DIR}/${prefix}.R2.paired.fastq.gz
    fastqc -t ${CPU} $file1 $file2
done < "$input"

# Set genome input files explicitly
GENOME_FA="${WD}/STAR_DB/input_files/genome.fa"
GENOME_GTF="${WD}/STAR_DB/input_files/annotation.gtf"

# Build index if not already present
if [ ! -f "${GENOME_DIR}/genomeParameters.txt" ]; then
    echo "Building STAR genome index..."
    mkdir -p "${GENOME_DIR}"

    STAR --runThreadN ${CPU} \
         --runMode genomeGenerate \
         --genomeDir ${GENOME_DIR} \
         --genomeFastaFiles ${GENOME_FA} \
         --sjdbGTFfile ${GENOME_GTF} \
         --sjdbOverhang 100
fi

# Align reads
ALIGN_OUT_DIR="${WD}/2-alignments"
mkdir -p "$ALIGN_OUT_DIR"

while IFS= read -r prefix; do
    echo "Aligning $prefix with STAR"
    STAR --runMode alignReads \
         --runThreadN ${CPU} \
         --genomeDir ${GENOME_DIR} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${ALIGN_OUT_DIR}/${prefix}. \
         --outSAMtype BAM SortedByCoordinate \
         --readFilesIn ${TRIM_DIR}/${prefix}.R1.paired.fastq.gz ${TRIM_DIR}/${prefix}.R2.paired.fastq.gz

    mv ${ALIGN_OUT_DIR}/${prefix}.Aligned.sortedByCoord.out.bam ${ALIGN_OUT_DIR}/${prefix}.sorted.bam
done < "$input"

input="./samples.prefix"
# Mark duplicates and generate BAM index

DEDUP_DIR="${WD}/3-deduplicated"
mkdir -p "$DEDUP_DIR"


while IFS= read -r prefix; do
    echo "Marking duplicates for $prefix"
    java -jar $PICARD_JAR MarkDuplicates \
        I=${ALIGN_OUT_DIR}/${prefix}.sorted.bam \
        O=${DEDUP_DIR}/${prefix}.sorted.dedup.bam \
        M=${DEDUP_DIR}/${prefix}.sorted.dedup.txt \
        READ_NAME_REGEX=null REMOVE_DUPLICATES=true

    samtools index ${DEDUP_DIR}/${prefix}.sorted.dedup.bam
done < "$input"


# Count reads with featureCounts
echo "Running featureCounts..."
files=()
while IFS= read -r prefix; do
    files+=(${DEDUP_DIR}/${prefix}.sorted.dedup.bam)
done < "$input"

COUNTS_DIR="${WD}/4-read_counts"
mkdir -p "$COUNTS_DIR"

featureCounts -p -a $GTF_ANNOTATION -T ${CPU} -s 0 -F GTF -t exon -g gene_id \
    -o ${COUNTS_DIR}/all_counts.txt -R BAM \
    --extraAttributes gene_name "${files[@]}"

module purge
