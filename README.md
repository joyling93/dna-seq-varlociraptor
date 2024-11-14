此流程仅用于WGS\WES\PANEL的DNA测序数据分析。

## 流程图示例
![流程图](./dag.svg "流程图示例")
## 流程环境
``conda activate /public/home/weiyifan/miniforge3/envs/sk8``
## 流程部署
``snakedeploy deploy-workflow https://github.com/joyling93/dna-seq-varlociraptor . --tag v1.1.0``
## 配置信息
### 必填
config.yaml  
samples.yaml  
units.yaml  
scenario.yaml  
### 选填
primers.tsv  
super_insteresting_genes.tsv  
## 流程运行
``snakemake -c30 --use-conda --cache``