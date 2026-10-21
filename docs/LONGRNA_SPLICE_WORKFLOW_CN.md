# NeoAg longRNA splice 工作流使用说明

本文说明如何使用 `scripts/run_longrna_splice_full.sh`，从长读长 RNA FASTQ 开始完成 transcript 重建、SQANTI3/TD2 注释、异常 splice junction 筛选、junction-spanning peptide 生成、MHC 预测和综合证据合并。

该流程用于研究级候选筛选和证据整理，不构成临床诊断或治疗建议。

## 1. 流程概览

完整流程依次执行：

1. 检查长读长 FASTQ 和必需输入。
2. 使用 IsoQuant 进行比对、转录本重建和表达统计。
3. 使用 SQANTI3 进行转录本结构分类和 junction 注释。
4. 通过 SQANTI3 的 `--include_ORF` 路径调用 TD2/PSAURON，生成 ORF 和蛋白质序列。
5. 排除短读长结果中已经存在的 junction。
6. 可选地排除正常背景队列中观察到的 junction。
7. 按 coverage、canonical、CDS、NMD 和 PSAURON 条件筛选候选。
8. 围绕 junction 生成 MHC-I peptide。
9. 运行已配置的 MHC presentation 预测器。
10. 合并来源、质量、长读长支持和 MHC 证据，并生成审计清单。

## 2. 软件和数据要求

建议运行环境：

- 64 位 Linux
- Bash 4 或更高版本
- Python 3.10 或更高版本
- IsoQuant 3.x
- SQANTI3 6.x 及其兼容 Python 环境
- SQANTI3 可调用的 TD2/PSAURON 环境
- NeoAg 及需要使用的 MHC 预测器
- 足够的 CPU、内存和磁盘空间

大型 Nanopore 数据建议至少使用 16 个 CPU 线程、64 GB 内存。`full` 保存模式可能保留大量 BAM、分块文件和缓存，应提前确认工作目录空间。

所有输入必须使用一致的参考版本，通常为 GRCh38。FASTA、GTF、minimap2 index、短读长 BAM 和 junction 坐标不能混用不同 genome build。

## 3. 必需输入

### 3.1 长读长 FASTQ 目录

`--input-dir` 指向一个样本的 FASTQ 文件目录，支持：

```text
*.fastq
*.fastq.gz
*.fq
*.fq.gz
```

目录中的这些文件会按文件名排序后共同交给 IsoQuant。不要在同一目录混放不同患者或不同样本的 FASTQ。

### 3.2 参考序列和注释

必需参数：

```text
--reference-fasta /refs/GRCh38.fa
--annotation-gtf /refs/gencode.gtf
```

如果已有与 FASTA 匹配的 minimap2 index，可提供：

```text
--reference-index /refs/GRCh38.mmi
```

### 3.3 HLA 文件

`--hla` 文件每行一个 HLA allele，例如：

```text
HLA-A*02:01
HLA-A*24:02
HLA-B*15:01
HLA-C*03:04
```

不要在该文件中加入患者姓名、密码、token 或其他敏感字段。

### 3.4 短读长 junction 比较文件

`--shortread-comparisons` 是必需 TSV，至少包含：

```text
junction_key
chr1:12345:12789:+
chr7:5521001:5521450:-
```

`junction_key` 使用：

```text
chromosome:左侧exon边界:右侧exon边界:strand
```

如果没有短读长验证数据，应提供只有表头的文件：

```text
junction_key
```

这种情况表示“未观察到可比较的短读长 junction”，不能表述为已经得到短读长阴性验证。

### 3.5 正常背景 junction 文件

`--normal-junctions` 可选，格式同样为包含 `junction_key` 的 TSV。候选若出现在该文件中，会在 peptide 生成前被排除。

如果不提供该参数，结果会明确记录：

```text
normal_background_status = UNASSESSED
```

不能将 `UNASSESSED` 解释为“正常样本中不存在”。

### 3.6 可选短读长支持

可向 IsoQuant 提供排序并索引的 matched RNA BAM：

```text
--illumina-bam /data/CASE001/rna.sorted.bam
```

可向 SQANTI3 提供 STAR junction coverage：

```text
--star-junctions /data/CASE001/SJ.out.tab
```

这两个输入用途不同，可以同时提供。

## 4. 运行前检查

克隆并进入仓库：

```bash
git clone https://github.com/open-neoag/neoag.git
cd neoag
```

确认入口脚本：

```bash
bash scripts/run_longrna_splice_full.sh --help
```

确认工具：

```bash
/path/to/isoquant --version
/path/to/sqanti-python --version
test -s /path/to/SQANTI3/sqanti3_qc.py
```

建议先使用 `--dry-run`。它会检查输入路径并打印 IsoQuant 命令，不启动重型计算。

## 5. 标准运行

`standard` 是默认保存模式，适合常规分析：

```bash
bash scripts/run_longrna_splice_full.sh \
  --sample-id CASE001 \
  --input-dir /data/CASE001/longrna_fastq \
  --workdir /work/CASE001/longrna \
  --reference-fasta /refs/GRCh38.fa \
  --reference-index /refs/GRCh38.mmi \
  --annotation-gtf /refs/gencode.gtf \
  --hla /data/CASE001/hla.txt \
  --shortread-comparisons /data/CASE001/shortread_junctions.tsv \
  --normal-junctions /refs/normal_junctions.tsv \
  --illumina-bam /data/CASE001/rna.sorted.bam \
  --star-junctions /data/CASE001/SJ.out.tab \
  --isoquant-bin /envs/longrna/bin/isoquant \
  --sqanti-home /tools/SQANTI3 \
  --sqanti-python /envs/sqanti3/bin/python \
  --threads 16 \
  --retention standard
```

不需要的可选参数可以删除。

## 6. 完整审计运行

需要保存尽可能完整的中间文件并计算 SHA-256 时，使用：

```bash
bash scripts/run_longrna_splice_full.sh \
  ...其他参数... \
  --retention full
```

`full` 模式会：

- 向 IsoQuant 添加 `--keep_tmp`。
- 保存命令执行顺序。
- 记录 Git HEAD 和工作区状态。
- 记录 Python、Conda 包和工具版本。
- 对登记的输入文件计算 SHA-256。
- 对工作目录中的输出和中间文件计算 SHA-256。
- 失败时保存日志和失败现场清单。

原始 FASTQ 和公共参考文件不会复制到工作目录；清单会记录其绝对路径、大小、修改时间和 SHA-256。

注意：`full` 模式的输入和输出哈希可能耗时较长，并会保留 IsoQuant 临时文件。第三方工具若在自身进程内部强制删除临时文件，工作流无法恢复已经被该工具删除的内容。

## 7. Dry-run

完整示例：

```bash
bash scripts/run_longrna_splice_full.sh \
  --sample-id CASE001 \
  --input-dir /data/CASE001/longrna_fastq \
  --workdir /work/CASE001/longrna \
  --reference-fasta /refs/GRCh38.fa \
  --annotation-gtf /refs/gencode.gtf \
  --hla /data/CASE001/hla.txt \
  --shortread-comparisons /data/CASE001/shortread_junctions.tsv \
  --isoquant-bin /envs/longrna/bin/isoquant \
  --sqanti-home /tools/SQANTI3 \
  --sqanti-python /envs/sqanti3/bin/python \
  --retention full \
  --dry-run
```

确认输出中的：

- `sample_id` 正确。
- FASTQ 文件数量符合预期。
- FASTA、GTF、index 和工作目录正确。
- `full` 模式的 IsoQuant 命令包含 `--keep_tmp`。
- 正常背景缺失时显示 `normal_background=UNASSESSED`。

## 8. 断点续跑

任务中断后，使用完全相同的输入和工作目录：

```bash
bash scripts/run_longrna_splice_full.sh \
  ...与首次运行相同的参数... \
  --resume
```

只有具有完成标志且关键输出仍存在的阶段才应复用。不要在同一个工作目录中更换：

- 样本 ID
- FASTQ 集合
- FASTA/GTF 或 genome build
- HLA 文件
- 短读长或正常背景定义
- peptide 长度和关键筛选参数

如需改变这些条件，应建立新的工作目录。

## 9. 主要输出

```text
WORKDIR/
├── isoquant/
├── sqanti3/
├── translation/
├── splice_candidates/
│   ├── longrna_unique_abnormal_junctions.tsv
│   ├── longrna_unique_junction_peptides.tsv
│   ├── longrna_unique_prediction_input.tsv
│   └── longrna_unique_candidate_summary.json
├── presentation/
├── integration/
│   ├── comprehensive_evidence.summary.json
│   └── verification.txt
├── provenance/
│   ├── commands.sh
│   ├── input_manifest.json
│   ├── output_manifest.json
│   └── failure_manifest.json
├── comprehensive_evidence.tsv
├── longrna_splice_status.tsv
├── README.workflow.txt
├── .isoquant_complete
└── .sqanti3_complete
```

`failure_manifest.json` 仅在失败时生成。

## 10. 候选筛选规则

默认策略包括：

- SQANTI3 abnormal structural category
- fusion 以外要求 `junction_category=novel`
- 不出现在短读长比较清单中
- 不出现在已提供的正常背景中
- unique junction coverage 至少为 3
- canonical junction
- complete CDS
- 非 predicted NMD
- PSAURON 至少为 0.8
- junction 位于 CDS 内
- 蛋白质序列只包含标准氨基酸

默认生成跨越 junction 的 `8,9,10,11mer`。可用以下参数修改：

```text
--peptide-lengths 8,9,10,11,12
```

支持范围为 8–14 aa。

## 11. 常见问题

### 找到多个 IsoQuant transcript-model GTF

每个工作目录原则上只运行一个样本。若自动检测到多个 GTF，应检查是否混入其他运行结果；通用 downstream 脚本也可以显式传入 `--isoform-gtf`。

### SQANTI3 在 isoform_hits 阶段报参数数量错误

仅对已经确认属于历史 API 冲突的 SQANTI3 版本使用：

```text
--patch-sqanti3
```

该参数会调用 `fix_sqanti3_isoform_hits_api.py`，并在修改前创建备份。不要把它作为所有 SQANTI3 版本的默认操作。

### 历史 SQANTI3 任务已经生成大部分结果但最后失败

可使用 `resume_sqanti3_outputs.py` 恢复对应历史失败类型。它不是常规新任务步骤。

### 候选数为零

零候选不一定是程序错误。首先查看：

```text
splice_candidates/longrna_unique_candidate_summary.json
```

重点检查 `excluded_reasons`、短读长重叠数、正常背景重叠数、coverage、CDS、NMD 和 PSAURON 条件。

### 综合证据合并返回缺失预测

检查 presentation 输出是否覆盖每个 `peptide + HLA allele` 组合，并确认没有跳过生产要求的预测器。缺失证据应保持显式，不应填充为阴性结果。

### 如何检查最终状态

```bash
bash scripts/verify_longrna_splice.sh /work/CASE001/longrna
```

成功时输出：

```text
workflow_status=PASS
```

## 12. Xiaoliang Site 1 历史兼容入口

历史样本 `260F400881011` 保留以下包装器：

```bash
bash scripts/run_longrna_downstream_260F400881011.sh
```

该脚本包含历史服务器路径，只适用于相同部署环境。新样本必须使用参数化入口：

```text
scripts/run_longrna_splice_full.sh
```

## 13. 引用与审计建议

正式报告至少保留：

- Git commit
- `provenance/input_manifest.json`
- `provenance/output_manifest.json`
- `provenance/commands.sh`
- `longrna_splice_status.tsv`
- 候选筛选摘要
- MHC predictor 原始结果和综合证据表
- 正常背景来源、版本和构建方式

如果结果用于跨批次比较，应固定参考版本、工具版本、HLA 输入、junction 坐标定义、正常背景和 peptide 长度。
