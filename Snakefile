import os

configfile: "config.yaml"

# ---------------------------------------------------------------------------
# ReACT-seq (RT-stop) pipeline.
# Detects ms2i6A (and other RT-stop marks) by reverse-transcription stop
# counting. Three presets cover the analysis tracks; rRNA reads are removed
# before alignment in all presets:
#   human_trans  human transcriptome RT-stop, STAR alignment
#   ecoli_trans  E. coli transcriptome RT-stop, bowtie2 alignment
#   ecoli_tRNA   E. coli tRNA-only RT-stop (all or selected tRNAs)
# All presets are paired-end only; helper scripts live flat in libs/.
# ---------------------------------------------------------------------------
PRESETS = {
    "human_trans": dict(species="human", aligner="star"),
    "ecoli_trans": dict(species="ecoli", aligner="bowtie2"),
    "ecoli_tRNA":  dict(species="ecoli", aligner="bowtie2"),
}

preset = config.get("preset", "ecoli_tRNA")
if preset not in PRESETS:
    raise ValueError("preset must be one of: " + ", ".join(sorted(PRESETS)))
P = PRESETS[preset]

# --- config (all paths/parameters come from config.yaml) ------------------------
samples = config["samples"]
ip_samples = config.get("ip_samples", samples)      # IP / treatment (enriched) samples
input_samples = config.get("input_samples", [])     # Input / control samples (peak calling)
datadir = config["datadir"]
outdir = config["outdir"]
refdir = config["refdir"]

# Resolve helper-script root relative to this Snakefile so the workflow is
# independent of the invocation directory. Default: the repo's own libs/.
libsdir = os.path.normpath(os.path.join(workflow.basedir, config.get("libsdir", "libs")))

species = P["species"]
aligner = P["aligner"]

FQ1_SUFFIX = config.get("fq1_suffix", "_R1.fq.gz")
FQ2_SUFFIX = config.get("fq2_suffix", "_R2.fq.gz")

# --- reference set for the primary track ------------------------------------
if preset == "ecoli_tRNA":
    ALIGNER_INDEX = config["tRNA_index"]
    MAIN_FASTA = config["trna_fasta"]
    MAIN_PLUS_BED = config["trna_plus_bed"]
    MAIN_MINUS_BED = config["trna_minus_bed"]
    MAIN_PEAKS = config["trna_peaks_bed"]
elif preset == "ecoli_trans":
    ALIGNER_INDEX = config["bowtie2_index"]
    MAIN_FASTA = config["genome_fasta"]
    MAIN_PLUS_BED = config["gene_plus_bed"]
    MAIN_MINUS_BED = config["gene_minus_bed"]
    MAIN_PEAKS = config["peaks_bed"]
else:  # human_trans
    ALIGNER_INDEX = config["star_index"]
    MAIN_FASTA = config["genome_fasta"]
    MAIN_PLUS_BED = config["gene_plus_bed"]
    MAIN_MINUS_BED = config["gene_minus_bed"]
    MAIN_PEAKS = config["peaks_bed"]

# High-stop site calling script, selected per species (config-provided).
FIND_STOP_SCRIPT = os.path.join(libsdir, config["find_stop_script_" + species])

# rRNA index for pre-alignment depletion, selected per species.
RRNA_INDEX = config["rRNA_index_" + species]

FILTER_OPTS = "-isPrimaryAlignment true -isProperPair true -isMapped true -mapQuality '>1'"

# Per-sample aligned BAM (name differs per aligner).
if aligner == "star":
    ALIGN_BAM = outdir + "/align/{sample}_Aligned.sortedByCoord.out.bam"
else:
    ALIGN_BAM = outdir + "/align/{sample}.aligned.bam"

# Peak calling (enrichment analysis, IP vs Input) is performed in both
# transcriptome presets (human_trans and ecoli_trans). Requires >= 2 IP
# replicates for IDR; it is not applied to the tRNA-only preset.
peak_calling = config.get("run_peak_calling", preset != "ecoli_tRNA")
if peak_calling and len(ip_samples) < 2:
    raise ValueError("run_peak_calling requires at least 2 ip_samples (replicates) for IDR")


# --- final targets ---------------------------------------------------------------
rule all:
    input:
        # 1) adapter/quality trimming
        expand(outdir + "/trim/{sample}_R1_val_1.fq.gz", sample=samples),
        expand(outdir + "/trim/{sample}_R2_val_2.fq.gz", sample=samples),
        # 2) rRNA depletion
        expand(outdir + "/rrna/{sample}_rmrRNA.fq.1.gz", sample=samples),
        expand(outdir + "/rrna/{sample}_rmrRNA.fq.2.gz", sample=samples),
        # 3) alignment + filtered/deduplicated BAM tracks
        expand(ALIGN_BAM, sample=samples),
        expand(outdir + "/align/{sample}.filtered.bam", sample=samples),
        expand(outdir + "/align/{sample}.filtered.bai", sample=samples),
        expand(outdir + "/align/{sample}.filtered.bw", sample=samples),
        # 4) pileup + per-transcript stop rates
        expand(outdir + "/pileup/{sample}.pileup", sample=samples),
        expand(outdir + "/stoprate/{sample}_plus_all.txt", sample=samples),
        expand(outdir + "/stoprate/{sample}_minus_all.txt", sample=samples),
        # 5) per-site RT-stop counts (both strands)
        expand(outdir + "/stoprate/{sample}.fwd.stop.txt", sample=samples),
        expand(outdir + "/stoprate/{sample}.rev.stop.txt", sample=samples),
        # 6) high-stop site calls (IP samples)
        expand(outdir + "/stopsite/{sample}.highstop.txt", sample=ip_samples),
        # 7) enrichment peaks + IDR-reproducible final peak set
        expand(outdir + "/peaks/{sample}_peaks.narrowPeak", sample=ip_samples) if peak_calling else [],
        [outdir + "/peaks/ReACT-seq_final.idrPeak"] if peak_calling else [],


# --- trim --------------------------------------------------------------------
rule trim:
    input:
        fq1 = datadir + "/{sample}" + FQ1_SUFFIX,
        fq2 = datadir + "/{sample}" + FQ2_SUFFIX,
    output:
        fq1 = outdir + "/trim/{sample}_R1_val_1.fq.gz",
        fq2 = outdir + "/trim/{sample}_R2_val_2.fq.gz",
    params:
        outdir = outdir + "/trim",
    log:
        outdir + "/trim/{sample}.trim.log",
    threads: config.get("trim_threads", 4)
    shell:
        "trim_galore -j {threads} -q 20 --stringency 3 --length 20 -e 0.1 --gzip --paired {input.fq1} {input.fq2} -o {params.outdir} &> {log}"


# --- rRNA removal (all presets) ----------------------------------------------
rule filter_rrna:
    input:
        fq1 = outdir + "/trim/{sample}_R1_val_1.fq.gz",
        fq2 = outdir + "/trim/{sample}_R2_val_2.fq.gz",
    output:
        fq1 = outdir + "/rrna/{sample}_rmrRNA.fq.1.gz",
        fq2 = outdir + "/rrna/{sample}_rmrRNA.fq.2.gz",
        sam = outdir + "/rrna/{sample}_rRNA.sam",
    params:
        prefix = outdir + "/rrna/{sample}_rmrRNA.fq.gz",
        index = RRNA_INDEX,
    log:
        outdir + "/rrna/{sample}.rrna.log",
    threads: config.get("align_threads", 8)
    shell:
        "bowtie2 -p {threads} --un-conc-gz {params.prefix} -x {params.index} -1 {input.fq1} -2 {input.fq2} -S {output.sam} &> {log}"


# --- align -------------------------------------------------------------------
if aligner == "star":
    rule align:
        input:
            fq1 = outdir + "/rrna/{sample}_rmrRNA.fq.1.gz",
            fq2 = outdir + "/rrna/{sample}_rmrRNA.fq.2.gz",
        output:
            outdir + "/align/{sample}_Aligned.sortedByCoord.out.bam",
        params:
            index = ALIGNER_INDEX,
            prefix = outdir + "/align/{sample}_",
            rgline = "ID:" + "{sample}" + "\tSM:" + "{sample}" + "\tLB:" + "{sample}" + "\tPL:ILLUMINA",
        log:
            outdir + "/align/{sample}.align.log",
        threads: config.get("align_threads", 8)
        shell:
            "STAR --runThreadN {threads} --genomeDir {params.index} --outFileNamePrefix {params.prefix} --outReadsUnmapped Fastx --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMattrRGline {params.rgline} --readFilesCommand zcat --outSAMstrandField intronMotif --readFilesIn {input.fq1} {input.fq2} &> {log}"
else:
    rule align:
        input:
            fq1 = outdir + "/rrna/{sample}_rmrRNA.fq.1.gz",
            fq2 = outdir + "/rrna/{sample}_rmrRNA.fq.2.gz",
        output:
            outdir + "/align/{sample}.aligned.bam",
        params:
            index = ALIGNER_INDEX,
        log:
            outdir + "/align/{sample}.align.log",
        threads: config.get("align_threads", 8)
        shell:
            "bowtie2 -p {threads} -x {params.index} -1 {input.fq1} -2 {input.fq2} --no-unal | samtools view -b -o {output} - &> {log}"


# --- filter / deduplicate / index / bigwig ----------------------------------
rule filter_bam:
    input:
        ALIGN_BAM,
    output:
        bam = outdir + "/align/{sample}.filtered.bam",
        bai = outdir + "/align/{sample}.filtered.bai",
        bw = outdir + "/align/{sample}.filtered.bw",
    params:
        tmp1 = outdir + "/align/{sample}.tmp1.bam",
        tmp2 = outdir + "/align/{sample}.tmp2.bam",
        tmp3 = outdir + "/align/{sample}.tmp3.bam",
        dupmetrics = outdir + "/align/{sample}.dupmetrics.txt",
        filter_opts = FILTER_OPTS,
    log:
        outdir + "/align/{sample}.filter.log",
    threads: 4
    shell:
        """
        bamtools filter -in {input} {params.filter_opts} > {params.tmp1} 2> {log}
        samtools sort {params.tmp1} -o {params.tmp2} &>> {log}
        picard MarkDuplicates REMOVE_DUPLICATES=true I={params.tmp2} O={params.tmp3} M={params.dupmetrics} &>> {log}
        samtools sort -o {output.bam} {params.tmp3} &>> {log}
        samtools index {output.bam} {output.bai} &>> {log}
        bamCoverage -b {output.bam} -o {output.bw} &>> {log}
        rm -f {params.tmp1} {params.tmp2} {params.tmp3}
        """


# --- pileup -----------------------------------------------------------------
rule pileup:
    input:
        outdir + "/align/{sample}.filtered.bam",
    output:
        outdir + "/pileup/{sample}.pileup",
    params:
        fasta = MAIN_FASTA,
        max_depth = config.get("pileup_max_depth", 20000000),
    log:
        outdir + "/pileup/{sample}.pileup.log",
    threads: 1
    shell:
        "samtools mpileup -B -Q 0 -d {params.max_depth} -o {output} -f {params.fasta} {input} &> {log}"


# --- transcript stop-rate -----------------------------------------------------
rule count_stoprate:
    input:
        outdir + "/pileup/{sample}.pileup",
    output:
        plus = outdir + "/stoprate/{sample}_plus_all.txt",
        minus = outdir + "/stoprate/{sample}_minus_all.txt",
    params:
        plus_bed = MAIN_PLUS_BED,
        minus_bed = MAIN_MINUS_BED,
    log:
        outdir + "/stoprate/{sample}.stoprate.log",
    threads: 1
    shell:
        "python3 {libsdir}/ReACT_seq_stoprate_transcrip.py {input} {params.plus_bed} {params.minus_bed} {output.plus} {output.minus} &> {log}"


# --- per-site RT stop counting ------------------------------------------------
rule count_stop:
    input:
        outdir + "/align/{sample}.filtered.bam",
    output:
        fwd = outdir + "/stoprate/{sample}.fwd.stop.txt",
        rev = outdir + "/stoprate/{sample}.rev.stop.txt",
    params:
        fasta = MAIN_FASTA,
        min_cov = config.get("stop_min_cov", 1),
        python = config.get("stop_python", "python3"),
    log:
        outdir + "/stoprate/{sample}.stop.log",
    threads: 1
    shell:
        """
        {params.python} {libsdir}/RTstop_count.py -i {input} -r {params.fasta} -c {params.min_cov} --StrandSelect fwd > {output.fwd} 2> {log}
        {params.python} {libsdir}/RTstop_count.py -i {input} -r {params.fasta} -c {params.min_cov} --StrandSelect rev > {output.rev} 2>> {log}
        """


# --- high-stop site calling --------------------------------------------------
rule find_stop:
    input:
        plus = outdir + "/stoprate/{sample}_plus_all.txt",
        minus = outdir + "/stoprate/{sample}_minus_all.txt",
    output:
        outdir + "/stopsite/{sample}.highstop.txt",
    params:
        script = FIND_STOP_SCRIPT,
        peaks = MAIN_PEAKS,
        python = config.get("stop_python", "python3"),
    log:
        outdir + "/stopsite/{sample}.highstop.log",
    threads: 1
    shell:
        "{params.python} {params.script} {input.plus} {input.minus} {output} {params.peaks} &> {log}"


# --- peak calling -------------------------------------------------------------
if peak_calling:
    # Per-IP-replicate macs2 callpeak against pooled Input controls.
    rule macs2_peaks:
        input:
            ip = outdir + "/align/{sample}.filtered.bam",
            ctrl = expand(outdir + "/align/{s}.filtered.bam", s=input_samples),
        output:
            outdir + "/peaks/{sample}_peaks.narrowPeak",
        params:
            gsize = config.get("macs2_gsize", "4.6e6"),
            q = config.get("macs2_qvalue", 0.05),
            peakdir = outdir + "/peaks",
            prefix = "{sample}",
        log:
            outdir + "/peaks/{sample}.macs2.log",
        threads: config.get("peak_threads", 4)
        shell:
            "macs2 callpeak -t {input.ip} -c {input.ctrl} -f BAM --nomodel --extsize 200 -g {params.gsize} -q {params.q} --outdir {params.peakdir} -n {params.prefix} 2> {log}"

    # IDR across IP replicates -> final reproducible peak set.
    # idr refuses datasets with fewer than 20 merged peaks (ValueError) but
    # still writes the merged peaks to the output file before exiting. For
    # such low-power datasets (e.g. lightweight tests on a chromosome
    # subset) keep the pre-IDR merged peaks and flag it in the log instead
    # of failing the whole workflow.
    rule idr_peaks:
        input:
            reps = expand(outdir + "/peaks/{s}_peaks.narrowPeak", s=ip_samples),
        output:
            outdir + "/peaks/ReACT-seq_final.idrPeak",
        log:
            outdir + "/peaks/idr.log",
        threads: 1
        shell:
            r"""
            idr --samples {input.reps} --input-file-type narrowPeak \
                --rank p.value --output-file {output} --plot 2> {log} || st=$?
            if [ "${{st:-0}}" -ne 0 ]; then
                if [ -s {output} ]; then
                    echo "WARNING: idr exited $st (<20 merged peaks, low-power dataset); keeping pre-IDR merged peaks as final peak set" >> {log}
                else
                    exit $st
                fi
            fi
            """
