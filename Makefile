
all: chr1.peaks.tab.gz chr2.peaks.tab.gz chr3.peaks.tab.gz chr4.peaks.tab.gz chr5.peaks.tab.gz chr6.peaks.tab.gz chr7.peaks.tab.gz chr8.peaks.tab.gz chr9.peaks.tab.gz chr10.peaks.tab.gz chr11.peaks.tab.gz chr12.peaks.tab.gz chr13.peaks.tab.gz chr14.peaks.tab.gz chr15.peaks.tab.gz chr16.peaks.tab.gz chr17.peaks.tab.gz chr18.peaks.tab.gz chr19.peaks.tab.gz chr20.peaks.tab.gz chr21.peaks.tab.gz chr22.peaks.tab.gz chrX.peaks.tab.gz

all_peaks.bed: raw/peaks_calling/*.bed
	cat raw/peaks_calling/*.bed | sort -k1,1 -k2,2n > all_peaks.bed

merged_peaks.bed: all_peaks.bed
	bedtools merge -i all_peaks.bed > merged_peaks.bed

variants_in_peaks.recode.vcf.gz: merged_peaks.bed /home/seb/masters/topmed/exome_modified_header_100000_chr_100s.vcf.gz
	vcftools --bed merged_peaks.bed --gzvcf /home/seb/masters/topmed/exome_modified_header_100000_chr_100s.vcf.gz --out variants_in_peaks --recode

homer.motifs:
	wget http://homer.ucsd.edu/homer/custom.motifs -O homer.motifs

HOCOMOCOv11_full_pwms_HUMAN_mono.txt:
	wget http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_pwms_HUMAN_mono.txt

HOCOMOCOv11_full_standard_thresholds_HUMAN_mono.txt:
	wget http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_standard_thresholds_HUMAN_mono.txt

hocomoco_thresholds.tab: HOCOMOCOv11_full_standard_thresholds_HUMAN_mono.txt
	cat $^ | cut -f 1,2 | grep -v "P-values" > $@

merged_peaks_1.bed: merged_peaks.bed
	awk '{print $$0 >> $$1".bed"}' merged_peaks.bed

variants_in_peaks_1.recode.vcf.gz: variants_in_peaks.recode.vcf.gz merged_peaks_1.bed
	vcftools --bed merged_peaks_1.bed --gzvcf variants_in_peaks.recode.vcf.gz --out variants_in_peaks --recode

ALL_BED=bed/Bcell-13.bed,bed/CD4-9.bed,bed/CD8-10.bed,bed/CLP-14.bed,bed/CMP-4.bed,bed/Erythro-15.bed,bed/GMP-5.bed,bed/HSC-1.bed,bed/LMPP-3.bed,bed/MCP.bed,bed/mDC.bed,bed/MEGA1.bed,bed/MEGA2.bed,bed/MEP-6.bed,bed/Mono-7.bed,bed/MPP-2.bed,bed/Nkcell-11.bed,bed/pDC.bed

# Not found: JFOS_HUMAN.H11MO.0.A,JFOSB_HUMAN.H11MO.0.A
MOTIFS=JUNB_HUMAN.H11MO.0.A,FOSL1_HUMAN.H11MO.0.A,FOSL2_HUMAN.H11MO.0.A,JDP2_HUMAN.H11MO.0.D,GATA1_HUMAN.H11MO.0.A,GATA2_HUMAN.H11MO.0.A,GATA3_HUMAN.H11MO.0.A,GATA4_HUMAN.H11MO.0.A,GATA5_HUMAN.H11MO.0.D,GATA6_HUMAN.H11MO.0.A,JUN_HUMAN.H11MO.0.A,JUND_HUMAN.H11MO.0.A,BATF_HUMAN.H11MO.0.A,ATF3_HUMAN.H11MO.0.A,BACH1_HUMAN.H11MO.0.A,BACH2_HUMAN.H11MO.0.A,NFE2_HUMAN.H11MO.0.A,CEBPA_HUMAN.H11MO.0.A,CEBPB_HUMAN.H11MO.0.A,CEBPD_HUMAN.H11MO.0.C,CEBPE_HUMAN.H11MO.0.A,CEBPG_HUMAN.H11MO.0.B,SPIB_HUMAN.H11MO.0.A,IRF8_HUMAN.H11MO.0.B,SPI1_HUMAN.H11MO.0.A,MESP1_HUMAN.H11MO.0.D,ID4_HUMAN.H11MO.0.D,HTF4_HUMAN.H11MO.0.A,ITF2_HUMAN.H11MO.0.C,STAT1_HUMAN.H11MO.0.A,STAT2_HUMAN.H11MO.0.A,SPIC_HUMAN.H11MO.0.D,CTCF_HUMAN.H11MO.0.A,IRF1_HUMAN.H11MO.0.A,DBP_HUMAN.H11MO.0.B,MAFK_HUMAN.H11MO.1.A,ATF4_HUMAN.H11MO.0.A,ASCL1_HUMAN.H11MO.0.A,ASCL2_HUMAN.H11MO.0.D,TFE2_HUMAN.H11MO.0.A,MYOD1_HUMAN.H11MO.0.A,EVI1_HUMAN.H11MO.0.B,IRF3_HUMAN.H11MO.0.B,ZEB1_HUMAN.H11MO.0.A,IRF9_HUMAN.H11MO.0.C,HEN1_HUMAN.H11MO.0.C,LYL1_HUMAN.H11MO.0.A

OPTIONS=--motifs $(MOTIFS) --bed-files $(ALL_BED) --ref-genome hg38.fa

COMMON_DEPS=hg38.fa HOCOMOCOv11_full_pwms_HUMAN_mono.txt hocomoco_thresholds.tab

chr1.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr1.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 1 --vcf raw/vcf/chr1.vcf.gz --output-file chr1.peaks.tab.gz

chr2.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr2.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 2 --vcf raw/vcf/chr2.vcf.gz --output-file chr2.peaks.tab.gz

chr3.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr3.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 3 --vcf raw/vcf/chr3.vcf.gz --output-file chr3.peaks.tab.gz

chr4.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr4.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 4 --vcf raw/vcf/chr4.vcf.gz --output-file chr4.peaks.tab.gz

chr5.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr5.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 5 --vcf raw/vcf/chr5.vcf.gz --output-file chr5.peaks.tab.gz

chr6.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr6.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 6 --vcf raw/vcf/chr6.vcf.gz --output-file chr6.peaks.tab.gz

chr7.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr7.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 7 --vcf raw/vcf/chr7.vcf.gz --output-file chr7.peaks.tab.gz

chr8.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr8.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 8 --vcf raw/vcf/chr8.vcf.gz --output-file chr8.peaks.tab.gz

chr9.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr9.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 9 --vcf raw/vcf/chr9.vcf.gz --output-file chr9.peaks.tab.gz

chr10.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr10.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 10 --vcf raw/vcf/chr10.vcf.gz --output-file chr10.peaks.tab.gz

chr11.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr11.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 11 --vcf raw/vcf/chr11.vcf.gz --output-file chr11.peaks.tab.gz

chr12.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr12.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 12 --vcf raw/vcf/chr12.vcf.gz --output-file chr12.peaks.tab.gz

chr13.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr13.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 13 --vcf raw/vcf/chr13.vcf.gz --output-file chr13.peaks.tab.gz

chr14.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr14.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 14 --vcf raw/vcf/chr14.vcf.gz --output-file chr14.peaks.tab.gz

chr15.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr15.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 15 --vcf raw/vcf/chr15.vcf.gz --output-file chr15.peaks.tab.gz

chr16.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr16.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 16 --vcf raw/vcf/chr16.vcf.gz --output-file chr16.peaks.tab.gz

chr17.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr17.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 17 --vcf raw/vcf/chr17.vcf.gz --output-file chr17.peaks.tab.gz

chr18.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr18.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 18 --vcf raw/vcf/chr18.vcf.gz --output-file chr18.peaks.tab.gz

chr19.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr19.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 19 --vcf raw/vcf/chr19.vcf.gz --output-file chr19.peaks.tab.gz

chr20.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr20.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 20 --vcf raw/vcf/chr20.vcf.gz --output-file chr20.peaks.tab.gz

chr21.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr21.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 21 --vcf raw/vcf/chr21.vcf.gz --output-file chr21.peaks.tab.gz

chr22.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chr22.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome 22 --vcf raw/vcf/chr22.vcf.gz --output-file chr22.peaks.tab.gz

chrX.peaks.tab.gz: $(COMMON_DEPS) raw/vcf/chrX.vcf.gz
	stack exec dnamotifs -- $(OPTIONS) --chromosome X --vcf raw/vcf/chrX.vcf.gz --output-file chrX.peaks.tab.gz