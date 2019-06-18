
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

ALL_BED=bed/Bcell-13.bed,bed/CD4-9.bed,bed/CD8-10.bed,bed/CLP-14.bed,bed/CMP-4.bed,bed/Erythro-15.bed,bed/GMP-5.bed,bed/HSC-1.bed,bed/LMPP-3.bed,bed/MCP.bed,bed/mDC.bed,bed/MEGA1.bed,bed/MEGA2.bed,bed/MEP-6.bed,bed/Mono-7.bed,bed/MPP-2.bed,bed/Nkcell-11.bed,bed/pDC.bed

chr1.peaks.tab.gz: hg38.fa raw/vcf/chr1.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 1 $(ALL_BED) $^ $@

chr2.peaks.tab.gz: hg38.fa raw/vcf/chr2.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 2 $(ALL_BED) $^ $@

chr3.peaks.tab.gz: hg38.fa raw/vcf/chr3.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 3 $(ALL_BED) $^ $@

chr4.peaks.tab.gz: hg38.fa raw/vcf/chr4.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 4 $(ALL_BED) $^ $@

chr5.peaks.tab.gz: hg38.fa raw/vcf/chr5.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 5 $(ALL_BED) $^ $@

chr6.peaks.tab.gz: hg38.fa raw/vcf/chr6.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 6 $(ALL_BED) $^ $@

chr7.peaks.tab.gz: hg38.fa raw/vcf/chr7.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 7 $(ALL_BED) $^ $@

chr8.peaks.tab.gz: hg38.fa raw/vcf/chr8.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 8 $(ALL_BED) $^ $@

chr9.peaks.tab.gz: hg38.fa raw/vcf/chr9.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 9 $(ALL_BED) $^ $@

chr10.peaks.tab.gz: hg38.fa raw/vcf/chr10.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 10 $(ALL_BED) $^ $@

chr11.peaks.tab.gz: hg38.fa raw/vcf/chr11.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 11 $(ALL_BED) $^ $@

chr12.peaks.tab.gz: hg38.fa raw/vcf/chr12.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 12 $(ALL_BED) $^ $@

chr13.peaks.tab.gz: hg38.fa raw/vcf/chr13.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 13 $(ALL_BED) $^ $@

chr14.peaks.tab.gz: hg38.fa raw/vcf/chr14.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 14 $(ALL_BED) $^ $@

chr15.peaks.tab.gz: hg38.fa raw/vcf/chr15.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 15 $(ALL_BED) $^ $@

chr16.peaks.tab.gz: hg38.fa raw/vcf/chr16.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 16 $(ALL_BED) $^ $@

chr17.peaks.tab.gz: hg38.fa raw/vcf/chr17.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 17 $(ALL_BED) $^ $@

chr18.peaks.tab.gz: hg38.fa raw/vcf/chr18.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 18 $(ALL_BED) $^ $@

chr19.peaks.tab.gz: hg38.fa raw/vcf/chr19.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 19 $(ALL_BED) $^ $@

chr20.peaks.tab.gz: hg38.fa raw/vcf/chr20.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 20 $(ALL_BED) $^ $@

chr21.peaks.tab.gz: hg38.fa raw/vcf/chr21.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 21 $(ALL_BED) $^ $@

chr22.peaks.tab.gz: hg38.fa raw/vcf/chr22.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- 22 $(ALL_BED) $^ $@

chrX.peaks.tab.gz: hg38.fa chrX.bed raw/vcf/chrX.vcf.gz HOCOMOCOv11_core_pwms_HUMAN_mono.txt hocomoco_thresholds.tab
	stack exec dnamotifs -- X $(ALL_BED) $^ $@