### RUN merge_fastq_lanes FIRST
### THEN RUN make merged_bam

.SECONDEXPANSION:
.SECONDARY:
.DELETE_ON_ERROR:

### SYSTEM TOOLS
SHELL = /bin/bash
LN = /bin/ln
MV = /bin/mv
CAT = /bin/cat
MKDIR = /bin/mkdir
ECHO = /bin/echo
CP = /bin/cp
CD = cd
SOURCE = source
AWK = /bin/awk
SORT = /bin/sort
UNIQ = /usr/bin/uniq
GZIP = /bin/gzip
ZCAT = /bin/zcat
GREP = /bin/grep
EGREP = /bin/egrep
SHUF = /usr/bin/shuf
SORT = /bin/sort
SED = /bin/sed
CUT = /bin/cut
PASTE = /usr/bin/paste

### USER TOOLS
BWA = /athena/abc/scratch/paz2005/bin/src/bwa-0.7.17/bwa #v0.7.17
SAMTOOLS = /athena/abc/scratch/paz2005/bin/src/samtools-1.8/samtools  #v 1.8
PYTHON = /home/paz2005/miniconda3/bin/python #v 3.6.5 
PYTHON2 = /home/paz2005/miniconda2/bin/python #v2.7.11
PYTHON3 = /home/paz2005/miniconda3/bin/python #v 3.6.5 
JAVA = /athena/abc/scratch/paz2005/bin/src/subread-1.6.2-Linux-x86_64/bin/jdk1.8.0_171/bin/java # v 1.8.0_171
PICARD = /athena/abc/scratch/paz2005/bin/src/picard-2.18.9/picard.jar #2.18.9
FASTQC = /athena/abc/scratch/paz2005/bin/src/FastQC/fastqc #  v0.11.7
RSCRIPT = /home/paz2005/miniconda3/bin/Rscript #v 3.4.3
R = /home/paz2005/miniconda3/bin/R #v 3.4.3
BEDTOOLS = /athena/abc/scratch/paz2005/bin/src/bedtools2/bin/bedtools #v2.27.1
MACS2 = /home/paz2005/miniconda2/bin/macs2 # 2.1.1.20160309
IDR = /home/paz2005/miniconda3/bin/idr #v 2.0.3
BAM_COVERAGE = /home/paz2005/miniconda3/bin/bamCoverage #v3.1.0
MULTI_BAM_SUMMARY = /home/paz2005/miniconda3/bin/multiBamSummary  #v3.1.0
PLOT_CORRELATION = /home/paz2005/miniconda3/bin/plotCorrelation  #v3.1.0
ANNOTATE_PEAKS = /home/paz2005/miniconda3/bin/annotatePeaks.pl #v4.9.1-pl5.22.0_5


### REFERENCES
ANNOTATION = /athena/abc/scratch/paz2005/references/GRCm38.p6/gencode.vM17.annotation.gtf
REFERENCE_DIR = /athena/abc/scratch/paz2005/references/GRCm38.p6/bwa
REFERENCE_PREFIX = GRCm38.primary_assembly.genome.fa
REFERENCE_FA = /athena/abc/scratch/paz2005/references/GRCm38.p6/GRCm38.primary_assembly.genome.fa
REFERENCE_INFO = /athena/abc/scratch/paz2005/references/GRCm38.p6/GRCm38.primary_assembly.genome.fa.fai
BLACKLIST_REGIONS = /athena/abc/scratch/paz2005/references/blacklist_regions/mm10.blacklist.merged.bed


### OPTIONS
CONDITIONS = N KO WT OT1 SV40
BWA_PIPELINE_OPTIONS = --threads 8
ORGANISM = mouse
MULTIMAPPING_NUMBER_CUTOFF = 4
NREADS=15000000
SHIFT_SIZE = -75
EXT_SIZE = 150
PVAL_THRES = 0.01
NEG_LOG10_QVAL_THRES = 2
IDR_THRESH = 0.1
BAM_COVERAGE_NRM_TO_1X = 2150570000  # mouse
EFFECTIVE_GENOME_SIZE = 1.87e9 # mouse
HOMER_ANNOTATION_CODE = mm10

### HElPER FUNCTIONS 
lc = $(subst A,a,$(subst B,b,$(subst C,c,$(subst D,d,$(subst E,e,$(subst F,f,$(subst G,g,$(subst H,h,$(subst I,i,$(subst J,j,$(subst K,k,$(subst L,l,$(subst M,m,$(subst N,n,$(subst O,o,$(subst P,p,$(subst Q,q,$(subst R,r,$(subst S,s,$(subst T,t,$(subst U,u,$(subst V,v,$(subst W,w,$(subst X,x,$(subst Y,y,$(subst Z,z,$1))))))))))))))))))))))))))

### DO NOT EDIT BELOW THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING
MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
CURRENT_DIR := $(patsubst %/,%,$(dir $(MKFILE_PATH)))
FASTQFILES := $(wildcard *_R1*.fastq.gz)
SAMPLES := $(sort $(foreach a,$(FASTQFILES),$(firstword $(subst _, ,$a))))
FLOWCELLS := $(sort $(join $(foreach a,$(FASTQFILES),$(firstword $(subst _, ,$a))), $(addprefix _, $(foreach a,$(FASTQFILES),$(word 2, $(subst _, ,$a))))))
BARCODES :=  $(sort $(join $(join $(foreach a,$(FASTQFILES),$(firstword $(subst _, ,$a))), $(addprefix _, $(foreach a,$(FASTQFILES),$(word 2, $(subst _, ,$a))))), $(addprefix _, $(foreach a,$(FASTQFILES),$(word 3, $(subst _, ,$a))))))
LANES := $(sort $(join $(join $(join $(foreach a,$(FASTQFILES),$(firstword $(subst _, ,$a))), $(addprefix _, $(foreach a,$(FASTQFILES),$(word 2, $(subst _, ,$a))))), $(addprefix _, $(foreach a,$(FASTQFILES),$(word 3, $(subst _, ,$a))))), $(addprefix _, $(foreach a,$(FASTQFILES),$(word 4, $(subst _, ,$a))))))
MERGED_BAMS := $(foreach a,$(SAMPLES), $(addsuffix $a.maxL.bam, $(addsuffix /,$a)) )
FASTQFILES := $(shell find *_R1_*.fastq.gz)
R2_FASTQFILES := $(wildcard *_R2_*.fastq.gz)  

### TARGETS
default: fastqc phase1 phase2 phase3 phase4
all: fastqc phase1 phase2 phase3 phase4
fastqc: $(patsubst %.fastq.gz,%_fastqc.zip,$(wildcard *_R1_*.fastq.gz))
merge_fastq_lanes: $(addsuffix _R1.fastq.gz, $(SAMPLES)) $(addsuffix _R2.fastq.gz, $(SAMPLES))   

phase1: final_bam big_wig
merged_bam: $(addsuffix _PF.maxL.bam, $(SAMPLES))
filtered_bam: $(addsuffix _PF.filt.srt.bam, $(SAMPLES))
mate_fixed_bam: $(addsuffix _PF.fixmate.bam, $(SAMPLES))
dupmarked_bam: $(addsuffix _PF.dupmark.bam, $(SAMPLES))
rmdup_bam: $(addsuffix _PF.rmdup.bam, $(SAMPLES))
final_nmsort_bam: $(addsuffix _PF.final.nmsrt.bam, $(SAMPLES))
rmdup_nmsort_bam: $(addsuffix _PF.rmdup.nmsrt.bam, $(SAMPLES))
rmdup_index: $(addsuffix _PF.rmdup.bam.bai, $(SAMPLES))
final_bam: $(addsuffix _PF.final.bam, $(SAMPLES))
big_wig: $(addsuffix _normTo1x.bw, $(SAMPLES))


phase2: tagalign shifted_tags
tagalign: $(addsuffix _PF.tagalign.gz, $(SAMPLES))
shifted_tags: $(addsuffix _PF.tn5.tagalign.gz, $(SAMPLES))

phase3: phase3_narrow
phase3_narrow: narrow_peaks filtered_narrow_peaks correlation_in_peaks annotated_narrow
narrow_peaks: $(addsuffix _PF.narrowPeak.gz, $(SAMPLES))
filtered_narrow_peaks: $(addsuffix _PF.narrowPeak.filt.gz, $(SAMPLES))
correlation_in_peaks: correlation_readsInPeaks_pearson.png
annotated_narrow: $(addsuffix _PF.narrowPeak.filt.homer.bed, $(SAMPLES))


### RUN FASTQC 
%_fastqc.zip: %.fastq.gz $$(subst R1,R2,%.fastq.gz)
	$(FASTQC) $^

%_fastqc.zip: %.fastq.gz
	$(FASTQC) $<

### MERGE FASTQ FILES BY SAMPLE NAME
find-fastq-files = $(sort $(filter $1_% , $(FASTQFILES)))
find-r2-fastq-files = $(sort $(filter $1_% , $(R2_FASTQFILES)))

define merge-fastq-files
$1_R1.fastq.gz: $(call find-fastq-files,$1)
	 $$(if $$(filter 1, $$(words $$^)),$$(LN) -fs $$^ $$@,$(CAT) $$^ >> $$@)
$1_R2.fastq.gz: $(call find-r2-fastq-files,$1)
	$$(if $$(filter 1, $$(words $$^)),$$(LN) -fs $$^ $$@,$(CAT) $$^ >> $$@)
endef

$(foreach s,$(SAMPLES),$(eval $(call merge-fastq-files,$s)))   

### PHASE 1
### ALIGN TO REFERENCE WITH BWA
find-fastq-files = $(sort $(filter $1_% , $(FASTQFILES)))

define align-bam-files
$1_PF.maxL.bam: $(call find-fastq-files,$1)
	$(BWA) aln -t 1 $(REFERENCE_DIR)/$(REFERENCE_PREFIX) $$(notdir $$(wildcard $1*_R1.fastq.gz)) > $$(basename $$(basename $$(wildcard $1*_R1.fastq.gz))).sai ; \
	$(BWA) aln -t 1 $(REFERENCE_DIR)/$(REFERENCE_PREFIX) $$(notdir $$(wildcard $1*_R2.fastq.gz )) > $$(basename $$(basename $$(wildcard $1*_R2.fastq.gz))).sai ; \
	$(BWA) sampe $(REFERENCE_DIR)/$(REFERENCE_PREFIX) $$(basename $$(basename $$(wildcard $1*_R1.fastq.gz))).sai $$(basename $$(basename $$(wildcard $1*_R2.fastq.gz))).sai $$(notdir $$(wildcard $1*_R1.fastq.gz)) $$(notdir $$(wildcard $1*_R2.fastq.gz)) > $1.sam ; \
	$(SAMTOOLS) view -bS $1.sam > $$@ ; \
	$(RM) $1.sam $$(basename $$(basename $$(wildcard $1*_R1.fastq.gz))).sai  $$(basename $$(basename $$(wildcard $1*_R2.fastq.gz))).sai 
endef

$(foreach s,$(SAMPLES),$(eval $(call align-bam-files,$s)))

### FILTER AND SORT BAM BY NAME
%_PF.filt.srt.bam: %_PF.maxL.bam
	$(SAMTOOLS) view -q 10 -F 524 -f 2 -u $< | $(SAMTOOLS) sort -n -O bam -T $<.tmp - > $@

### FIX MATE COORDINATES THEN FILTER
%_PF.fixmate.bam: %_PF.filt.srt.bam
	$(SAMTOOLS) fixmate -r $< - | $(SAMTOOLS) view -F 1804 -f 2 -u - | $(SAMTOOLS) sort - > $@

### MARK DUPLICATES
%_PF.dupmark.bam: %_PF.fixmate.bam
	$(JAVA) -Xmx4g -jar $(PICARD) MarkDuplicates VALIDATION_STRINGENCY=SILENT INPUT=$< OUTPUT=$@ METRICS_FILE=$(basename $@).metrics REMOVE_DUPLICATES=false

### CALCULATE PERCENT DUPLICATES
%_PF.dupmark.log: %_PF.dupmark.bam
	$(SED) -n '/LIBRARY/,/^$$/p' $(basename $<).metrics | $(CUT) -f7,9 | $(EGREP) [0-9] | $(AWK) '{OFS="\t"; print "duplicated pairs:",$$1,$$2*100"%"}' > $@

### REMOVE DUPLICATES (final)
%_PF.rmdup.bam: %_PF.dupmark.bam
	$(SAMTOOLS) view -F 1804 -f 2 -b $< > $@

### INDEX BAM
%_PF.rmdup.bam.bai: %_PF.rmdup.bam
	$(SAMTOOLS) index $<

### NAME SORT FINAL BAM
%_PF.rmdup.nmsrt.bam: %_PF.rmdup.bam
	$(SAMTOOLS) sort -n -O bam -T $<.tmp $< > $@

#### CREATE FINAL BAM WITH CHRM AND SCAFFOLDS REMOVED
%_PF.final.bam: %_PF.rmdup.bam
	$(SAMTOOLS) view -h $< | egrep 'chr[0-9XY]+|^@HD|^@PG|XT:A:U'  | $(SAMTOOLS) view -bS - | $(SAMTOOLS) sort -n -O bam - > $@.tmp ; $(SAMTOOLS) fixmate -r $@.tmp - | $(SAMTOOLS) view -F 1804 -f 2 -u - | $(SAMTOOLS) sort - > $@ && rm $@.tmp

### NAME SORT FINAL BAM
%_PF.final.nmsrt.bam: %_PF.final.bam
	$(SAMTOOLS) sort -n -O bam -T $<.tmp $< > $@

### INDEX BAM
%_PF.final.bam.bai: %_PF.final.bam
	$(SAMTOOLS) index $<

### CREATE BIG WIG FILE
%_normTo1x.bw: %_PF.final.bam %_PF.final.bam.bai
	$(BAM_COVERAGE) -b $< -o $@ -bs 10 --normalizeUsing RPGC --effectiveGenomeSize $(BAM_COVERAGE_NRM_TO_1X) --blackListFileName $(BLACKLIST_REGIONS) --ignoreForNormalization chrX chrY --ignoreDuplicates --minFragmentLength 40 -p 1

### PHASE 2
### CONVERT BAM TO TAG ALIGN
%_PF.tagalign.gz: %_PF.final.bam
	$(BEDTOOLS) bamtobed -i $< | $(AWK) 'BEGIN{OFS="\t"}{$$4="N";$$5="1000";print $0}' | $(GZIP) -c > $@

### TN5 SHIFT TAG ALIGNS
%_PF.tn5.tagalign.gz: %_PF.tagalign.gz
	$(ZCAT) $< | $(AWK) -F $$'\t' 'BEGIN {OFS = FS}{ if ($$6 == "+") {$$2 = $$2 + 4} else if ($$6 == "-") {$$3 = $$3 - 5} print $$0}' | $(GZIP) -c > $@


### CALL NARROW PEAKS
%_PF_peaks.narrowPeak: %_PF.tn5.tagalign.gz
	$(MACS2) callpeak -t $< -f BED -n $(basename $(basename $(basename $^))) -g $(EFFECTIVE_GENOME_SIZE) -p $(PVAL_THRES) --nomodel --shift $(SHIFT_SIZE) --extsize $(EXT_SIZE) --SPMR --keep-dup all 
# -B --call-summits

%_PF.narrowPeak.gz: %_PF_peaks.narrowPeak
	$(SORT) -k 8gr,8gr $< | $(AWK) 'BEGIN{OFS="\t"}{$$4="Peak_"NR ; print $$0}' | $(GZIP) -c > $@

%_PF.narrowPeak.filt.gz: %_PF.narrowPeak.gz
	$(BEDTOOLS) intersect -v -a $< -b $(BLACKLIST_REGIONS) | $(AWK) 'BEGIN{OFS="\t"} {if ($$5>1000) $$5=1000; print $$0}' | $(GREP) -P 'chr[0-9XY]+(?!_)' |  $(AWK) 'BEGIN{OFS="\t"} $$9 >= $(NEG_LOG10_QVAL_THRES)' | $(BEDTOOLS) sort -i stdin | $(GZIP) -c > $@

%_PF.narrowPeak.filt.homer.bed: %_PF.narrowPeak.filt.gz
	$(ANNOTATE_PEAKS) <(zcat $<) $(HOMER_ANNOTATION_CODE) -gtf $(ANNOTATION) -annStats $(basename $@).annotation_stats.txt > $@

### PLOT CORRELATION BASED ON READS IN PEAKS
peaks_union.bed: $(addsuffix _PF.narrowPeak.filt.gz, $(SAMPLES))
	$(ZCAT) $^ | $(CUT) -f 1-3 | $(SORT) | $(UNIQ) | $(BEDTOOLS) sort -i stdin | $(BEDTOOLS) merge -i stdin > $@
	
readsInPeaks.npz: $(addsuffix _PF.final.bam, $(SAMPLES)) peaks_union.bed
	$(MULTI_BAM_SUMMARY) BED-file -b $(filter-out peaks_union.bed,$^) --BED peaks_union.bed -out $@ --outRawCounts $(basename $@).txt --ignoreDuplicates --blackListFileName $(BLACKLIST_REGIONS) --minFragmentLength 40 --maxFragmentLength 170

correlation_readsInPeaks_pearson.png: readsInPeaks.npz
	$(PLOT_CORRELATION) --corData $< -o $@ -p heatmap -c pearson --skipZeros --plotNumbers --removeOutliers --zMin .5


