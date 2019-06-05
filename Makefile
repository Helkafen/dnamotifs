
all: chr1.peaks.tab chr2.peaks.tab chr3.peaks.tab chr4.peaks.tab chr5.peaks.tab chr6.peaks.tab chr7.peaks.tab chr8.peaks.tab chr9.peaks.tab chr10.peaks.tab chr11.peaks.tab chr12.peaks.tab chr13.peaks.tab chr14.peaks.tab chr15.peaks.tab chr16.peaks.tab chr17.peaks.tab chr18.peaks.tab chr19.peaks.tab chr20.peaks.tab chr21.peaks.tab chr22.peaks.tab chrX.peaks.tab

all_peaks.bed: raw/peaks_calling/*.bed
	cat raw/peaks_calling/*.bed | sort -k1,1 -k2,2n > all_peaks.bed

merged_peaks.bed: all_peaks.bed
	bedtools merge -i all_peaks.bed > merged_peaks.bed

variants_in_peaks.recode.vcf.gz: merged_peaks.bed /home/seb/masters/topmed/exome_modified_header_100000_chr_100s.vcf.gz
	vcftools --bed merged_peaks.bed --gzvcf /home/seb/masters/topmed/exome_modified_header_100000_chr_100s.vcf.gz --out variants_in_peaks --recode

homer.motifs:
	wget http://homer.ucsd.edu/homer/custom.motifs -O homer.motifs

HOCOMOCOv11_core_pwms_HUMAN_mono.txt:
	wget http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_pwms_HUMAN_mono.txt

HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt:
	wget http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt

hocomoco_thresholds.tab: HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt
	cat $^ | cut -f 1,2 | grep -v "P-values" > $@

merged_peaks_1.bed: merged_peaks.bed
	awk '{print $$0 >> $$1".bed"}' merged_peaks.bed

variants_in_peaks_1.recode.vcf.gz: variants_in_peaks.recode.vcf.gz merged_peaks_1.bed
	vcftools --bed merged_peaks_1.bed --gzvcf variants_in_peaks.recode.vcf.gz --out variants_in_peaks --recode

chr%.peaks.tab: chr%.bed raw/vcf/chr%.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 1 hg38.fa $^ $@
