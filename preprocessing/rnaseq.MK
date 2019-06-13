.SECONDEXPANSION:
.SECONDARY:
.DELETE_ON_ERROR:

### SYSTEM TOOLS
LN = /bin/ln
MV = /bin/mv
CAT = /bin/cat
MKDIR = /bin/mkdir
ECHO = /bin/echo
CP = /bin/cp
CD = cd

### USER TOOLS
STAR = /athena/abc/scratch/paz2005/bin/src/STAR-2.6.0c/bin/Linux_x86_64_static/STAR #2.6.0c
TRIM_GALORE = /athena/abc/scratch/paz2005/bin/src/TrimGalore-0.5.0/trim_galore # v 0.5.0
SAMTOOLS = /athena/abc/scratch/paz2005/bin/src/samtools-1.8/samtools  #v 1.8
FEATURE_COUNTS = /athena/abc/scratch/paz2005/bin/src/subread-1.6.2-Linux-x86_64/bin/featureCounts #v1.6.2
JAVA = /athena/abc/scratch/paz2005/bin/src/subread-1.6.2-Linux-x86_64/bin/jdk1.8.0_171/bin/java # v 1.8.0_171
FASTQC = /athena/abc/scratch/paz2005/bin/src/FastQC/fastqc #  v0.11.7
QORTS = /athena/abc/scratch/paz2005/bin/src/QoRTs/QoRTs.jar # v1.3.0
RSCRIPT = /home/paz2005/miniconda3/bin/Rscript #v 3.4.3
R = /home/paz2005/miniconda3/bin/R #v 3.4.3

### REFERENCES
REFERENCE = /athena/abc/scratch/paz2005/references/GRCm38.p6
REFERENCE_FA =/athena/abc/scratch/paz2005/references/GRCm38.p6/GRCm38.primary_assembly.genome.fa
ANNOTATION = /athena/abc/scratch/paz2005/references/GRCm38.p6/gencode.vM17.annotation.gtf

### PARAMETERS
STAR_OPTIONS = --runThreadN 1 \
--runMode alignReads \
--genomeDir $(REFERENCE) \
--readFilesCommand zcat \
--outSAMstrandField intronMotif  \
--outFilterIntronMotifs RemoveNoncanonicalUnannotated  \
--outFilterType BySJout \
--outReadsUnmapped None \
--outSAMtype BAM Unsorted \
--chimOutType SeparateSAMold \
--limitSjdbInsertNsj 10000000 


### DO NOT EDIT BELOW THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING
FASTQFILES :=  $(wildcard *_R1*.fastq.gz)
SAMPLES := $(sort $(foreach a,$(FASTQFILES),$(firstword $(subst _, ,$a))))

all: fastqc merged_bam gene_counts qc_plots
fastqc: $(patsubst %.fastq.gz,%_fastqc.zip,$(wildcard *_R1_???.fastq.gz))
first_pass: $(FASTQFILES:.fastq.gz=.1p.Aligned.out.bam) $(FASTQFILES:.fastq.gz=.1p.SJ.out.tab)
sorted_bam: $(FASTQFILES:.fastq.gz=.sorted.bam)
merged_bam: $(addsuffix .merged.bam, $(SAMPLES))
gene_counts: $(addsuffix .gene.counts, $(SAMPLES))
qc_data: $(addprefix ./qc/, $(addsuffix /QC.summary.txt, $(SAMPLES)))


### FASTQC 
%_fastqc.zip: %.fastq.gz $$(subst R1,R2,%.fastq.gz)
	$(FASTQC) $^

### FIRST PASS ALIGNMENT
%.1p.Aligned.out.bam %.1p.SJ.out.tab: %.fastq.gz $$(subst R1,R2,%.fastq.gz)
	$(STAR) --genomeLoad LoadAndKeep --readFilesIn $^ --outFileNamePrefix $*.1p. $(STAR_OPTIONS)

### SORT BAM FILES
%.sorted.bam: %.1p.Aligned.out.bam
	$(SAMTOOLS) sort -n -O bam -T $<.tmp $< > $@

### MERGE BAM FILES BY SAMPLE NAME
find-bam-files = $(sort $(filter $1_% , $(FASTQFILES:.fastq.gz=.sorted.bam)))

define merge-bam-files
$1.merged.bam: $(call find-bam-files,$1)
	$$(if $$(filter 1, $$(words $$^)),$$(LN) -fs $$^ $$@,$(SAMTOOLS) merge -f $$@ $$^)
endef

$(foreach s,$(SAMPLES),$(eval $(call merge-bam-files,$s)))

### COUNT READS PER GENE
%.gene.counts : %.merged.bam
	file_arr=(./$**R2*.fastq.gz) && if [ -f $${file_arr[0]} ] ; then $(FEATURE_COUNTS) -p -g "gene_name" -a $(ANNOTATION) -o $@  $< ; else $(FEATURE_COUNTS) -g "gene_name" -a $(ANNOTATION) -o $@  $< ; fi

### QC
./qc/%/QC.summary.txt: %.merged.bam
	file_arr=(./$*_*R2*.fastq.gz) && $(MKDIR) -p ./qc/$* && raw=$$(echo $$(zcat $*_*R1*.fastq.gz | wc -l)/4 | bc) && if [ -f $${file_arr[0]}  ] ; then $(JAVA) -Xmx120G -jar $(QORTS) QC --seqReadCt $$(echo $$raw) $< $(ANNOTATION) ./qc/$* ; else $(JAVA) -Xmx120G -jar $(QORTS) QC --singleEnded --seqReadCt $$(echo $$raw) $< $(ANNOTATION) ./qc/$* ; fi
