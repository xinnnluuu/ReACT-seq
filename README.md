# ReACT-seq

> [中文说明 (Chinese)](README.zh-CN.md)

**Paper**: The Transcriptome-Wide Mapping of 2-Methylthio-N6-isopentenyladenosine (ms²i⁶A) at Single-Base Resolution by ReACT-seq
**Journal**: Journal of the American Chemical Society (2023). DOI: 10.1021/jacs.2c13894 (PMID: 36820840)

ReACT-seq maps ms²i⁶A across the E. coli (and human) transcriptome at
single-base resolution. Detection relies on **RT-stop**: the reverse
transcriptase terminates ahead of the modified site, and the modification is
localised by counting stop signals. This is an enrichment-style experiment
(IP vs Input); the workflow provides RT-stop site analysis and, in the
transcriptome presets, macs2 + IDR enrichment peak calling. rRNA reads are
removed before alignment in every preset.

## Directory layout
```
ReACT-seq/
├── README.md
├── Snakefile + config.yaml    # main workflow (Snakemake, 3 presets, PE only)
├── test_config.yaml           # chr21-subset test configuration
└── libs/                      # RT-stop core scripts
    ├── ReACT_seq_stoprate_transcrip.py   # per-transcript/gene RT-stop rate
    ├── RTstop_count.py                   # per-site stop counting
    ├── ecoli_RT_findstopsite.py          # E. coli high-stop site calling
    └── find_highstop_hg38.py             # human high-stop site calling
```

## Key presets (config.yaml -> preset)
| preset | description |
|---|---|
| ecoli_tRNA | E. coli **tRNA-only** RT-stop: rRNA removal -> align to tRNA index -> stop rate -> high-stop sites. tRNA reference can be all tRNAs or a selected subset |
| ecoli_trans | E. coli transcriptome RT-stop: rRNA removal -> transcriptome alignment -> stop rate -> high-stop sites + peak calling |
| human_trans | human transcriptome RT-stop: rRNA removal -> STAR alignment -> stop rate -> high-stop sites + peak calling |

All presets are **paired-end** and remove rRNA first; the workflow does not
include expression quantification (htseq-count / RPKM).

## Analysis steps
1. **trim** — Trim Galore (paired-end);
2. **rRNA removal** — bowtie2 depletion of rRNA reads (all presets);
3. **align** — bowtie2 (ecoli_tRNA: tRNA index; ecoli_trans: transcriptome index; human_trans: STAR);
4. **filter** — bamtools (primary/proper/mapped/mapQ>1) + Picard dedup + bigWig;
5. **pileup** — samtools mpileup at every site;
6. **stoprate** — ReACT_seq_stoprate_transcrip.py: RT-stop rate per transcript / tRNA;
7. **count_stop** — RTstop_count.py: per-site stop counts per strand;
8. **find_highstop** — ecoli_RT_findstopsite.py / find_highstop_hg38.py: high-stop site calling (ms²i⁶A candidates);
9. (transcriptome presets) **peak calling** — macs2 callpeak (each IP replicate vs pooled Input) -> IDR across replicates -> `ReACT-seq_final.idrPeak`.

## Usage
```bash
# 1) edit config.yaml: preset / samples / reference paths
# 2) dry-run, then run
snakemake -s Snakefile --configfile config.yaml -n
snakemake -s Snakefile --configfile config.yaml -j 8
```

## Reference file formats (config.yaml)
| key | format |
|---|---|
| star_index | STAR index directory (from `STAR --runMode genomeGenerate`) |
| bowtie2_index / tRNA_index / rRNA_index_* | bowtie2 index prefix (the part before `.1.bt2`) |
| genome_fasta / trna_fasta | FASTA (tRNA may be all or a subset, consistent with the BEDs) |
| gene_plus_bed / gene_minus_bed | BED3+name: `chr start end gene_id` (plus/minus strand, 0-based) |
| trna_plus_bed / trna_minus_bed | BED3+name: `chr start end trna_id` (plus/minus strand) |
| peaks_bed / trna_peaks_bed | BED3 candidate regions: `chr start end` (restricts high-stop calling) |

## Software
snakemake, trim_galore, cutadapt, bowtie2, STAR (human_trans), bamtools,
samtools, picard, deeptools (bamCoverage), macs2 + IDR (transcriptome
presets), python3, R (exomePeak2, optional).

## Data
GEO: GSE221930 (E. coli, treated vs control, replicates; see the GEO page for
SRA / BioProject accessions).

## Test note (test_config.yaml)
- Test data: ERR1122128/29 (IP) and ERR1122130 (Input) from GSE221930 on an
  hg38 chr21-subset reference — for **process validation only**, not
  genome-wide results.
- Validated end-to-end (2026-10): all steps (trim / rRNA removal / STAR
  alignment / filtering & dedup / pileup / stop-rate / per-site stops /
  high-stop sites / macs2 peaks / IDR) complete on the trio above.
- The IDR step requires >= 20 peaks post-merge; on the chr21 subset fewer
  peaks overlap, so `idr_peaks` keeps idr's pre-IDR merged peaks and records
  a `WARNING: ... low-power dataset` line in `peaks/idr.log` instead of
  failing (idr still writes the merged peaks before its size check).
  Full-genome datasets are unaffected.
- References for the test are rebuilt with `bash test_ref_prep.sh` (chr21
  fasta + strand-split gene BEDs + STAR/bowtie2 indexes + candidate peak
  regions from the workflow's own IP narrowPeaks).
