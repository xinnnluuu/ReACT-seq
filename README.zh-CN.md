# ReACT-seq

> [English](README.md)

**论文**: The Transcriptome-Wide Mapping of 2-Methylthio-N6-isopentenyladenosine (ms²i⁶A) at Single-Base Resolution by ReACT-seq
**期刊**: Journal of the American Chemical Society (2023). DOI: 10.1021/jacs.2c13894 (PMID: 36820840)

ReACT-seq 对大肠杆菌转录组的 2-甲基硫代-N6-异戊烯腺苷 (ms²i⁶A) 进行单碱基分辨全转录组图谱绘制。
检测原理为 **RT-stop (逆转录终止)**: 逆转录酶在修饰位点前终止, 通过统计 stop 信号定位修饰。
方法为富集型实验 (IP vs Input), 本流程提供 RT-stop 位点分析, 并在转录组预设中进行 macs2 + IDR 富集峰分析。
所有预设均先去 rRNA 再进行比对。

## 目录结构
```
ReACT-seq/
├── README.md
├── Snakefile + config.yaml    # 主流程 (Snakemake, 3 预设, 仅双端)
├── test_config.yaml           # chr21 子集测试配置
└── libs/                      # RT-stop 核心脚本
    ├── ReACT_seq_stoprate_transcrip.py   # 转录本/基因 RT-stop 率
    ├── RTstop_count.py                   # 逐位点 stop 计数
    ├── ecoli_RT_findstopsite.py          # E. coli 高 stop 位点筛选
    └── find_highstop_hg38.py             # human 高 stop 位点筛选
```

## 关键预设 (config.yaml -> preset)
| preset | 说明 |
|---|---|
| ecoli_tRNA | 大肠杆菌 **tRNA-only** RT-stop: 去除 rRNA 后比对到 tRNA 索引 → stoprate → 高 stop 位点。tRNA 参考可为全部 tRNA 或指定子集 |
| ecoli_trans | 大肠杆菌转录组 RT-stop: 去除 rRNA 后比对转录组 → stoprate → 高 stop 位点 + peak calling |
| human_trans | 人转录组 RT-stop: 去除 rRNA 后 STAR 比对 → stoprate → 高 stop 位点 + peak calling |

所有预设均为**双端**测序, 均先去 rRNA; 流程不含基因表达定量 (htseq-count/RPKM) 步骤。

## 分析步骤
1. **trim** — Trim Galore (双端);
2. **rRNA removal** — bowtie2 去除 rRNA 读段 (所有预设);
3. **align** — bowtie2 比对 (ecoli_tRNA: tRNA 索引; ecoli_trans: 转录组索引; human_trans: STAR);
4. **filter** — bamtools 过滤 primary/proper/mapped/mapQ>1 + Picard 去重 + bigWig;
5. **pileup** — samtools mpileup 全位点覆盖;
6. **stoprate** — ReACT_seq_stoprate_transcrip.py 计算每转录本/每 tRNA 的 RT-stop 率;
7. **count_stop** — RTstop_count.py 按链统计 stop 计数;
8. **find_highstop** — ecoli_RT_findstopsite.py / find_highstop_hg38.py 筛选高 stop 位点 (ms²i⁶A 候选);
9. (转录组预设) **peak calling** — macs2 callpeak (每个 IP 重复 vs pooled Input) → IDR 合并重复 → `ReACT-seq_final.idrPeak`。

## 用法
```bash
# 1) 修改 config.yaml: preset / samples / 参考文件路径
# 2) 预演与实跑
snakemake -s Snakefile --configfile config.yaml -n
snakemake -s Snakefile --configfile config.yaml -j 8
```

## 参考文件格式 (config.yaml)
| 配置项 | 格式 |
|---|---|
| star_index | STAR 索引目录 (STAR --runMode genomeGenerate 生成) |
| bowtie2_index / tRNA_index / rRNA_index_* | bowtie2 索引前缀 (`.1.bt2` 之前的部分) |
| genome_fasta / trna_fasta | FASTA (tRNA 可为全部或指定子集, 与 BED 保持一致) |
| gene_plus_bed / gene_minus_bed | BED3+name: `chr start end gene_id` (正/负链, 0-based) |
| trna_plus_bed / trna_minus_bed | BED3+name: `chr start end trna_id` (正/负链) |
| peaks_bed / trna_peaks_bed | BED3 候选区域: `chr start end` (限制高 stop 位点调用范围) |

## 所需软件
snakemake, trim_galore, cutadapt, bowtie2, STAR (human_trans), bamtools, samtools,
picard, deeptools (bamCoverage), macs2 + IDR (转录组预设), python3, R (exomePeak2, 可选)。

## 数据
GEO: GSE221930 (E. coli, 处理 vs 对照, 多重复; SRA/BioProject 见 GEO 页面)。

