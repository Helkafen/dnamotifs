
all_peaks.bed: raw/*.bed
	cat raw/*.bed | sort -k1,1 -k2,2n > all_peaks.bed

merged_peaks.bed: all_peaks.bed
	bedtools merge -i all_peaks.bed > merged_peaks.bed

variants_in_peaks.recode.vcf.gz: merged_peaks.bed /home/seb/masters/topmed/exome_modified_header_100000_chr_100s.vcf.gz
	vcftools --bed merged_peaks.bed --gzvcf /home/seb/masters/topmed/exome_modified_header_100000_chr_100s.vcf.gz --out variants_in_peaks --recode

homer.motifs:
	wget http://homer.ucsd.edu/homer/custom.motifs -O homer.motifs

merged_peaks_1.bed: merged_peaks.bed
	awk '{print $0 >> $1".bed"}' merged_peaks.bed

variants_in_peaks_1.recode.vcf.gz: variants_in_peaks.recode.vcf.gz merged_peaks_1.bed
	vcftools --bed merged_peaks_1.bed --gzvcf variants_in_peaks.recode.vcf.gz --out variants_in_peaks --recode

matches_1.tab: homer.motifs merged_peaks_1.bed variants_in_peaks_1.recode.vcf.gz
	stack exec dnamotifs -- 1 homer.motifs hg38.fa merged_peaks_1.bed variants_in_peaks_1.recode.vcf.gz matches_1.tab
