#!/bin/bash

#download data
#50 meters
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001326/SRR001326.fastq.gz
#1 meter
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001325/SRR001325.fastq.gz
#32 meters
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001323/SRR001323.fastq.gz
#16 meters
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001322/SRR001322.fastq.gz
#uncompress data
gunzip -v *.fastq.gz

#assembly - preparatory phase
velveth velvet_work 31 -fmtAuto  *.fastq
#assembly
velvetg velvet_work -min_contig_lgth 50
#change sequence names
cat velvet_work/contigs.fa | sed -E 's/(>NODE_[0-9]+)_.+/\1/g' > assembly.fa

#translation of the assembly in AA
translate_seq assembly.fa assembly.aa.fa

#download offline data
download_data -p -m EMAIL

#download profile for the three taxa
download_profiles -o bacteria_profiles -i 200795 201174 203682 -m EMAIL -k mg_data/kegg.pickle -t mg_data/taxonomy.pickle
#and archaea
download_profiles -o archaea_order_profiles -r order -l archaea -m EMAIL -k mg_data/kegg.pickle -t mg_data/taxonomy.pickle
download_profiles -o archaea_phylum_profiles -r phylum -l archaea -m EMAIL -k mg_data/kegg.pickle -t mg_data/taxonomy.pickle

#Profile alignments and build
for file in bacteria_profiles/*.fa; do
        clustalo -i $file -o bacteria_profiles/`basename $file .fa`.afa;
        hmmbuild bacteria_profiles/`basename $file .afa`.hmm $file;
done
for file in archaea_order_profiles/*.fa; do
        clustalo -i $file -o archaea_order_profiles/`basename $file .fa`.afa;
        hmmbuild archaea_order_profiles/`basename $file .afa`.hmm $file;
done
for file in archaea_phylum_profiles/*.fa; do
        clustalo -i $file -o archaea_phylum_profiles/`basename $file .fa`.afa;
        hmmbuild archaea_phylum_profiles/`basename $file .afa`.hmm $file;
done
#Make one file for all profiles
cat profile_files/*.hmm archaea_order_profiles/*.hmm archaea_phylum_profiles/*.hmm > profile_files.hmm

#HMMER command line
hmmsearch -o /dev/null --domtbl hmmer_dom-table.txt profile_files.hmm assembly.aa.fa
#convert into annotations
hmmer2gff -o assembly.gff -k mg_data/kegg.pickle -n assembly.fa assembly.aa.fa hmmer_dom-table.txt
#Filter annotations
filter_gff -l -i assembly.gff -o assembly.filt.gff

#extract aa sequences
#python -c "import mgkit;mgkit.logger.config_log(); import mgkit.io.gff; mgkit.io.gff.extract_aa_seqs('assembly.filt.gff', 'assembly.aa.pred.fa')"
#blastp
#blastp -query assembly.aa.pred.fa -db nr -out assembly.aa.pred.blast.tab -outfmt 7 -num_threas 8

#bowtie2
bowtie2-build assembly.fa assembly.fa
for file in *.fastq; do
	BASENAME=`basename $file .fastq`
	bowtie2 -N 1 -x assembly.fa -U $file \
	--rg-id $BASENAME --rg PL:454 --rg PU:454 \
	--rg SM:$BASENAME| samtools view -Sb - > $BASENAME.bam;
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
#.dict index;
#if there's an error execute: execstack -c picard-tools/libIntelDeflater.so
java -jar picard-tools/CreateSequenceDictionary.jar R=assembly.fa O=assembly.dict

#add coverage data
export SAMPLES=$(for file in *.bam; do echo -a `basename $file -sort.bam`,$file ;done)
add_coverage_to_gff -v $SAMPLES assembly.filt.gff assembly.filt.cov.gff
unset SAMPLES

#SNP calling using GATK
# #GATK
# for file in *-sort.bam; do

# 	  java -jar GATK/GenomeAnalysisTK.jar -R assembly.fa \
# 	  -T UnifiedGenotyper -I $file \
# 	  -o `basename $file .bam`.vcf -dcov 50;

# done

#SNP calling using samtools 
for file in *.bam; do 
	samtools mpileup -Iuf assembly.fa \
	$file |bcftools view -vcg - > `basename $file`.vcf;
done

#merge vcf
export SAMPLES=$(for file in *-sort.bam.vcf; do echo -V:`basename $file -sort.bam.vcf` $file ;done)
java -Xmx10g -jar GATK/GenomeAnalysisTK.jar \
	  -R assembly.fa -T CombineVariants -o assembly.vcf \
	  -genotypeMergeOptions UNIQUIFY \
	  $SAMPLES
unset SAMPLES

#convert to gtf
python -c "import mgkit;mgkit.logger.config_log(); import mgkit.io.gff; mgkit.io.gff.convert_gff_to_gtf('assembly.filt.gff', 'assembly.filt.gtf')"

#Download SNPdat
https://snpdat.googlecode.com/files/SNPdat_v1.0.5.pl

#SNPdat
for file in *SR*.vcf; do 
	perl SNPdat_v1.0.5.pl -i $file -g assembly.filt.gtf \
	-f assembly.fa -o `basename $file .vcf`.snpdat;
done

#snp_parser
export SAMPLES=$(for file in *.bam; do echo `basename $file -sort.bam`;done)
export SNPDAT_FILES=*.snpdat
snp_parser -v -b -g assembly.filt.cov.gff \
	-p assembly.vcf \
 	-s $SNPDAT_FILES \
 	-m $SAMPLES \
 	-t mg_data/taxonomy.pickle
unset SAMPLES
unset SNPDAT_FILES
