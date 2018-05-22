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
velveth velvet_work 31 -fmtAuto *.fastq
#assembly
velvetg velvet_work -min_contig_lgth 50
#change sequence names
cat velvet_work/contigs.fa | sed -E 's/(>NODE_[0-9]+)_.+/\1/g' > assembly.fa
#alternative
#fasta-utils uid cat velvet_work/contigs.fa assembly.fa
#remove velvet working directory
rm -R velvet_work

#To use the LCA option and other analysis we need a taxonomy file
download-taxonomy.sh

#Uniprot SwissProt DB
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
#Uncompress it
gunzip uniprot_sprot.fasta.gz

########
#Gene prediction

###BLAST
#index Uniprot
makeblastdb -dbtype prot -in uniprot_sprot.fasta
#Run blastx
blastx -query assembly.fa -db uniprot_sprot.fasta -out \
	assembly.uniprot.tab -outfmt 6

###RAPSearch
#Index
# prerapsearch -d uniprot_sprot.fasta -n uniprot_sprot
#Run
# rapsearch -q assembly.fa -d uniprot_sprot -o assembly.uniprot.tab
#rename .m8 file to assembly.uniprot.tab and delete .aln
#rm assembly.uniprot.tab.aln
#mv assembly.uniprot.tab.m8 assembly.uniprot.tab

########
#Converts gene prediction into GFF annotations
blast2gff uniprot -b 40 -db UNIPROT-SP -dbq 10 assembly.uniprot.tab \
	assembly.uniprot.gff
filter-gff overlap assembly.uniprot.gff assembly.uniprot-filt.gff
#rename the new filtered file
mv assembly.uniprot-filt.gff assembly.uniprot.gff

########
#Taxonomic refinement - requires NCBI nt DB installed and indexed
export NCBINT_DIR=ncbi-db
if [ -d "$NCBINT_DIR" ]; then
	echo "Taxonomic refinement";
	#Extract annotations sequences
	get-gff-info sequence -f assembly.fa assembly.uniprot.gff \
        assembly.uniprot.frag.fasta
	#Use blastn to match against NCBI NT
	blastn -query assembly.uniprot.frag.fasta -db ncbi-db/nt -out \
		assembly.uniprot.frag.tab -outfmt 6

	#Download necessary data
	download-ncbi-taxa.sh
	gunzip -c ncbi-nucl-taxa.gz | taxon-utils to_hdf -n nt

	# Get the LCA for the sequences
	blast2gff blastdb -i 3 -r assembly.uniprot.frag.tab | \
        filter-gff sequence -t -a length -f quantile -l 0.95 -c gt | \
        filter-gff sequence -t -a identity -f quantile -l 0.95 -c gt | \
        add-gff-info addtaxa -f taxa-table.hf5:nt | \
        taxon-utils lca -b 40 -t taxonomy.pickle -s -p - lca.tab

	# Add the LCA info to the GFF
	add-gff-info addtaxa -v -t lca.tab -a seq_id -db NCBI-NT \
        assembly.uniprot.gff assembly.uniprot-taxa.gff

	#rename the file to continue the script
	mv assembly.uniprot-taxa.gff assembly.uniprot.gff
fi
unset NCBINT_DIR

########
#Finalise information from Gene Prediction
#Adds remaining taxonomy, EC numbers, KO and eggNOG mappings
add-gff-info uniprot --buffer 500 -t -e -ec -ko \
	assembly.uniprot.gff assembly.uniprot-final.gff
#Rename the GFF
mv assembly.uniprot-final.gff assembly.uniprot.gff

########
#Alignments
bowtie2-build assembly.fa assembly.fa
for file in *.fastq; do
	BASENAME=`basename $file .fastq`;
	bowtie2 -N 1 -x assembly.fasta -U $file \
	--very-sensitive-local \
	--rg-id $BASENAME --rg PL:454 --rg PU:454 \
	--rg SM:$BASENAME | samtools view -Sb - > $BASENAME.bam;
done
#sort and index BAM files with samtools
for file in *.bam; do
	samtools sort $file `basename $file .bam`-sort;
	#removes the unsorted file, it's not needed
	mv `basename $file .bam`-sort.bam $file
	samtools index $file;
done

########
#Add coverage and expected changes to GFF file
#export SAMPLES=$(for file in *.bam; do echo -a `basename $file .bam`,$file ;done)
#Coverage info
#add-gff-info coverage $SAMPLES assembly.uniprot.gff | add-gff-info \
#	exp_syn -r assembly.fa > assembly.uniprot-update.gff
#rename to continue the script
#mv assembly.uniprot-update.gff assembly.uniprot.gff
#unset SAMPLES
#Faster
for x in *.bam; do samtools depth -aa $x > `basename $x .bam`.depth; done
add-gff-info cov_samtools $(for file in *.depth; do echo -s `basename file .depth` -d $file ;done) assembly.uniprot.gff assembly.uniprot-update.gff
mv assembly.uniprot-update.gff assembly.uniprot.gff

########
#SNP calling using samtools
bcftools mpileup -f assembly.fa -Ou *.bam | bcftools call -v -O v --ploidy 1 -A -o assembly.vcf


#Index fasta with Picard tools - GATK requires it
#java -jar picard-tools/picard.jar CreateSequenceDictionary \
#	R=assembly.fa O=assembly.dict

#merge vcf
#export SAMPLES=$(for file in *.bam.vcf; do echo -V:`basename $file .bam.vcf` $file ;done)
# java -Xmx10g -jar GATK/GenomeAnalysisTK.jar \
# 	  -R assembly.fa -T CombineVariants -o assembly.vcf \
# 	  -genotypeMergeOptions UNIQUIFY \
# 	  $SAMPLES
# unset SAMPLES

########
#snp_parser
export SAMPLES=$(for file in *.bam; do echo -m `basename $file .bam`;done)
snp_parser -s -v -g assembly.uniprot.gff -p assembly.vcf -a assembly.fa $SAMPLES
unset SAMPLES

########
#Count reads
# Using HTSeq
for file in *.bam; do
	htseq-count -f bam -r pos -s no -t CDS -i uid -a 8 \
	-m intersection-nonempty $file assembly.uniprot.gff \
	> `basename $file .bam`-counts.txt
done

# Adding counts to GFF
add-gff-info counts `for x in *.bam; do echo -s $(basename $x .bam); done` \
	`for x in *-counts.txt; do echo -c $x; done` assembly.uniprot.gff tmp.gff

mv tmp.gff assembly.uniprot.gff
# Using featureCounts
#featureCounts -a assembly.uniprot.gff -g uid -O  -t CDS -o counts-featureCounts.txt *.bam
#add-gff-info counts `for x in *.bam; do echo -s $(basename $x .bam); done` -c counts-featureCounts.txt -e assembly.uniprot.gff tmp.gff
#mv tmp.gff assembly.uniprot.gff

########
#eggNOG mappings
wget http://eggnog.embl.de/version_3.0/data/downloads/COG.members.txt.gz
wget http://eggnog.embl.de/version_3.0/data/downloads/NOG.members.txt.gz
wget http://eggnog.embl.de/version_3.0/data/downloads/COG.funccat.txt.gz
wget http://eggnog.embl.de/version_3.0/data/downloads/NOG.funccat.txt.gz

########
#EC names
wget ftp://ftp.expasy.org/databases/enzyme/enzclass.txt
