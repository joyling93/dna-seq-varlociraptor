# 一般设置
要配置此工作流，请根据你的需求修改“config/config.yaml”，并遵循文件中提供的说明。

# 样本表
将样本添加到“config/samples.tsv”中。对于每个样本，必须定义“sample_name”、“alias”、“platform”、“datatype”、“calling”和“group”列。
  - 同一“group”中的样本可以通过它们的“alias”在联合[call 突变场景](#calling-scenario)中被引用。
  - “alias”表示样本在其组中的名称。它们旨在作为样本类型的一些抽象描述，用于[call 突变场景](#calling-scenario)，因此应在不同组中一致使用。一个典型的例子是“tumor”和“normal”别名的组合。
  - “platform”列需要包含使用的测序平台（“CAPILLARY”、“LS454”、“ILLUMINA”、“SOLID”、“HELICOS”、“IONTORRENT”、“ONT”、“PACBIO”之一）。
  - 在使用默认场景时需要纯度列。如果未知，可以设置为“1.0”。
  - 在“samples.tsv”样本表中，相同的“sample_name”条目可以多次使用，只有重复行之间“group”列的值不同。这样，你可以在不同的组中使用相同的样本 call 突变，例如，当你没有匹配的正常样本用于肿瘤call 突变时，可以使用一组正常样本。
  - “datatype”列指定每个样本对应的数据类型。可以是“rna”或“dna”。
  - “calling”列设置要执行的分析类型。可以是“fusions”、“variants”或两者（用逗号分隔）。
如果要估计样本的突变负荷，则必须在一个额外的列“mutational_burden_events”中指定来自call 突变场景的要使用的“事件”（见下文）。该列中的多个事件必须用逗号分隔。
缺失值可以通过空列或写入“NA”来指定。行可以用“#”注释掉。

# 单元表
对于每个样本，将一个或多个测序单元（运行、分组或重复）添加到单元表“config/units.tsv”中。
  - 每个单元有一个“unit_name”。这可以是一个连续的编号，也可以是实际的运行、分组或重复 ID。
  - 每个单元有一个“sample_name”，它将单元与它来自的生物样本相关联。此信息用于在读取映射和重复标记之前合并样本的所有单元。
  - 对于每个单元，你需要指定以下列中的一个：
    - 仅对于单端读取指定“fq1”。这可以指向系统上的任何 FASTQ 文件。
    - 对于成对末端读取指定“fq1”和“fq2”。这些可以指向系统上的任何 FASTQ 文件。
    - 仅指定“sra”：指定一个 SRA（序列读取档案）登录号（例如以 ERR 或 SRR 开头）。管道将自动从 SRA 下载相应的成对末端读取。
    - 如果本地文件（“fq1”、“fq2”）和 SRA 登录号（“sra”）都可用，则使用本地文件。
  - 在“adapters”列中定义适配器，将[cutadapt 参数](https://cutadapt.readthedocs.org)放在引号中（例如“-a ACGCGATCG -A GCTAGCGTACT”）。
缺失值可以通过空列或写入“NA”来指定。行可以用“#”注释掉。

# call 突变场景
Varlociraptor 支持对任意场景的变异进行集成的不确定性感知call 突变和过滤。这些被定义为所谓的场景，通过[变异call 突变语法](https://varlociraptor.github.io/docs/calling#generic-variant-calling)。
  - [scenario.yaml 模板参考](https://varlociraptor.github.io/varlociraptor-scenarios/landing/)
  - [格式说明](https://varlociraptor.github.io/docs/calling/)

# 引物修剪
当启用引物修剪时，引物必须直接在“config.yaml”中或在单独的 tsv 文件中定义。
当所有样本来自同一引物组时，首选在配置文件中直接定义引物。
对于不同的pannel，引物必须在单独的 tsv 文件中按设置。
对于每个pannel，需要设置以下列：“panel”、“fa1”和“fa2”（可选）。
此外，对于每个样本，相应的pannel必须在“samples.tsv”中定义（“panel”列）。
仅对于单引物修剪，必须定义配置（或相应的 tsv 文件）中的第一个条目。

# 注释 UMIs
为了注释 UMIs，必须在“sample.tsv”中设置两个额外的列：
  - “umi_read”：可以是以下选项之一：
    - 如果 UMIs 是读取 1 的一部分，则为“fq1”。
    - 如果 UMIs 是读取 2 的一部分，则为“fq2”。
    - 如果成对末端读取中都有 UMIs，则为“both”。
    - 指向一个额外的 fastq 文件的路径，该文件仅包含 fq1 和 fq2 中每个片段的 UMI（具有相同的读取名称）。
  - “umi_read_structure”：定义每个 UMI 记录中 UMI 位置的读取结构（见 https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures）。如果两个读取都包含一个 UMI，则在两者之间用空格指定一个读取结构（例如，“8M142T 8M142T”）。如果设置了仅包含 UMI 序列的单独 fastq 文件，则读取结构必须为“+M”。
UMI 记录的读取名称必须与样本 fastqs 的相应记录匹配。