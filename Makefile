
all: chr1.peaks.tab.gz chr2.peaks.tab.gz chr3.peaks.tab.gz chr4.peaks.tab.gz chr5.peaks.tab.gz chr6.peaks.tab.gz chr7.peaks.tab.gz chr8.peaks.tab.gz chr9.peaks.tab.gz chr10.peaks.tab.gz chr11.peaks.tab.gz chr12.peaks.tab.gz chr13.peaks.tab.gz chr14.peaks.tab.gz chr15.peaks.tab.gz chr16.peaks.tab.gz chr17.peaks.tab.gz chr18.peaks.tab.gz chr19.peaks.tab.gz chr20.peaks.tab.gz chr21.peaks.tab.gz chr22.peaks.tab.gz chrX.peaks.tab.gz

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

chr1.peaks.tab.gz: hg38.fa chr1.bed raw/vcf/chr1.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 1 $^ $@

chr2.peaks.tab.gz: hg38.fa chr2.bed raw/vcf/chr2.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 2 $^ $@

chr3.peaks.tab.gz: hg38.fa chr3.bed raw/vcf/chr3.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 3 $^ $@

chr4.peaks.tab.gz: hg38.fa chr4.bed raw/vcf/chr4.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 4 $^ $@

chr5.peaks.tab.gz: hg38.fa chr5.bed raw/vcf/chr5.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 5 $^ $@

chr6.peaks.tab.gz: hg38.fa chr6.bed raw/vcf/chr6.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 6 $^ $@

chr7.peaks.tab.gz: hg38.fa chr7.bed raw/vcf/chr7.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 7 $^ $@

chr8.peaks.tab.gz: hg38.fa chr8.bed raw/vcf/chr8.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 8 $^ $@

chr9.peaks.tab.gz: hg38.fa chr9.bed raw/vcf/chr9.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 9 $^ $@

chr10.peaks.tab.gz: hg38.fa chr10.bed raw/vcf/chr10.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 10 $^ $@

chr11.peaks.tab.gz: hg38.fa chr11.bed raw/vcf/chr11.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 11 $^ $@

chr12.peaks.tab.gz: hg38.fa chr12.bed raw/vcf/chr12.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 12 $^ $@

chr13.peaks.tab.gz: hg38.fa chr13.bed raw/vcf/chr13.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 13 $^ $@

chr14.peaks.tab.gz: hg38.fa chr14.bed raw/vcf/chr14.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 14 $^ $@

chr15.peaks.tab.gz: hg38.fa chr15.bed raw/vcf/chr15.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 15 $^ $@

chr16.peaks.tab.gz: hg38.fa chr16.bed raw/vcf/chr16.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 16 $^ $@

chr17.peaks.tab.gz: hg38.fa chr17.bed raw/vcf/chr17.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 17 $^ $@

chr18.peaks.tab.gz: hg38.fa chr18.bed raw/vcf/chr18.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 18 $^ $@

chr19.peaks.tab.gz: hg38.fa chr19.bed raw/vcf/chr19.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 19 $^ $@

chr20.peaks.tab.gz: hg38.fa chr20.bed raw/vcf/chr20.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 20 $^ $@

chr21.peaks.tab.gz: hg38.fa chr21.bed raw/vcf/chr21.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 21 $^ $@

chr22.peaks.tab.gz: hg38.fa chr22.bed raw/vcf/chr22.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 22 $^ $@

chrX.peaks.tab.gz: hg38.fa chrX.bed raw/vcf/chrX.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- X $^ $@