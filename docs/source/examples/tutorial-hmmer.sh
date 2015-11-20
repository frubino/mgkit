#!/bin/bash

# Some tools require a contact email
export EMAIL=your@email

# Requirements
#
# software 	- version tested
# mgkit 	- 0.2.1
# wget 		- any
# velvet 	- 1.2.10
# bowtie2 	- 2.2.6 / 2.1.0
# khmer 	- 2.0
# hmmer 	- 3.1b2 / 3.1b1
# clustalo 	- 1.2.1

# Requirements specific for Mac OS X (El Capitan, 10.11): wget
# brew install wget velvet bowtie2 samtools pyenv-virtualenv hmmer clustal-omega

# Memory ~6 GB for khmer normalisation
# ~9 GB for patitioning
# ~3 GB for velveth
# ~1.5 GB for velvetg
# ~3.5 GB for megahit

# Using data from the following project:
# http://www.ebi.ac.uk/ena/data/view/PRJEB6461
# Samples:
# I (Influent)
# B (Buffering)
# SA (Secondary aeration)
# PA (Primary aeration)
# SD (Sludge digestion)
mkdir seqs
cd seqs
#download data
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/B_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/B_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/I_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/I_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/PA_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/PA_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/SA_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/SA_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/SD_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/SD_R2.fastq.gz

# Normalisation of the reads
# Khmer when using paired reads works best when are interleaved
# Moreover, it requires the older Casava type of header
# The sequences are trimmed by 20 bp and fed to the khmer
# script interleave-reads.py
interleave-reads.py --gzip -o all-interleaved.fq.gz \
<(
python - <<END
import HTSeq, sys, glob
files = glob.glob('*R1.fastq.gz')
for fname in files:
    for record in HTSeq.FastqReader(fname):
        record = record[:-20] # cut 20 bp
        # HTSeq adds '[part]' to the header, the next line cut it
        record.name = record.name[:-6]
        record.write_to_fastq_file(sys.stdout)
END
) \
<(
python - <<END
import HTSeq, sys, glob
files = glob.glob('*R2.fastq.gz')
for fname in files:
    for record in HTSeq.FastqReader(fname):
        record = record[:-20]
        record.name = record.name[:-6]
        record.write_to_fastq_file(sys.stdout)
END
)
# Normalisation takes out almost 7% of the reads
normalize-by-median.py -k 24 -p -o normalised.fq.gz --gzip -M 6e9 all-interleaved.fq.gz
# rm all-interleaved.fq.gz

megahit --presets meta --verbose -t 2 --min-contig-len 50 --12 normalised.fq.gz -o megahit-out

cd ..

# Get the nitrogen metabolism KOs
export GENES=`python - <<END
from mgkit import kegg
k = kegg.KeggClientRest()
print " ".join(k.link_ids('ko', 'ko00910')['ko00910'])
END`

#download offline data (only taxonomy)
download_data -p -x -m $EMAIL

download_profiles -o profiles -ko $GENES -r order -l bacteria -m $EMAIL -t mg_data/taxonomy.pickle

unset GENES

#Profile alignments and build
for file in profiles/*.fa; do
	echo $file `basename $file .fa`
    clustalo -i $file -o profiles/`basename $file .fa`.afa ;
    hmmbuild profiles/`basename $file .afa`.hmm profiles/`basename $file .fa`.afa;
done
#Make one file for all profiles
cat profiles/*.hmm > profiles.hmm

#translation of the assembly in AA
translate_seq assembly.fa assembly.aa.fa

#HMMER command line
hmmsearch -o /dev/null --domtbl hmmer_dom-table.txt profiles.hmm assembly.aa.fa
#convert into annotations
hmmer2gff -d -o assembly.gff assembly.aa.fa hmmer_dom-table.txt
#Filter annotations
filter-gff overlap -s 100 assembly.gff assembly.filt.gff

#bowtie2
bowtie2-build assembly.fa assembly.fa
for file in R1.fastq.gz; do
	BASENAME=`basename $file _R1.fastq.gz`
	bowtie2 -N 1 -x assembly.fa -1 $file -2 $file_R2.fastq.gz \
	--rg-id $BASENAME --rg PL:Illumina --rg PU:Illumina-MiSeq \
	--rg SM:$BASENAME | samtools view -Sb - > $BASENAME.bam;
done

#samtools
for file in *.bam; do
	samtools sort $file `basename $file .bam`-sort;
	samtools index `basename $file .bam`-sort.bam;
	#removes the unsorted file, it's not needed
	rm $file;
done

#index for the assembly (required by GATK and samtools)
#.fai index
samtools faidx assembly.fa

#add coverage data
export SAMPLES=$(for file in *.bam; do echo -a `basename $file -sort.bam`,$file ;done)
add_coverage_to_gff -v $SAMPLES assembly.filt.gff assembly.filt.cov.gff

#SNP calling using samtools
for file in *.bam; do
	samtools mpileup -Iuvf assembly.fa $file > `basename $file`.vcf;
done

samtools mpileup -t DP,DV -s -q 10 -Iuvf assembly.fa *.bam

#fasta .dict file, needed for GATK
wget https://github.com/broadinstitute/picard/releases/download/1.140/picard-tools-1.140.zip
unzip picard-tools-1.140.zip
java -jar picard-tools-1.140/picard.jar CreateSequenceDictionary \
	R=assembly.fa O=assembly.dict

#merge vcf
export SAMPLES=$(for file in *-sort.bam.vcf; do echo -V:`basename $file -sort.bam.vcf` $file ;done)
java -Xmx1g -jar GATK/GenomeAnalysisTK.jar \
	  -R assembly.fa -T CombineVariants -o assembly.vcf \
	  -genotypeMergeOptions UNIQUIFY \
	  $SAMPLES
unset SAMPLES

#snp_parser
export SAMPLES=$(for file in *.bam; do echo -m `basename $file .bam`;done)
snp_parser -v -g assembly.uniprot.gff \
	-p assembly.vcf \
	-a assembly.fasta \
 	$SAMPLES
unset SAMPLES
