import pandas as pd
antb = snakemake.input[0]
database=pd.read_csv('/public/home/xiezhuoming/xzm/ref/vcf/human/CGC/vcf_database_gene.tsv',header=0,sep="\t")
annotation=pd.read_csv(antb,header=0,sep="\t")

driver=pd.merge(annotation.iloc[:,0],database,left_on = 'symbol',right_on='Gene Symbol', how='inner')
driver.drop_duplicates().to_csv(snakemake.output[0],sep="\t")