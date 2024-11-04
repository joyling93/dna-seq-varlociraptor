if config["purity_ploidy"]["activate"]:
    rule pp_prep:
        input:
            "results/recal/{sample}.bam"
        output:
            hdo = "results/purity_ploidy/{sample}.hdo.sam",
            hdn = "results/purity_ploidy/{sample}.hdn.sam",
            bam = "results/purity_ploidy/{sample}.chr.bam",
            idx = "results/purity_ploidy/{sample}.chr.bam.bai",
        shell:
            """
            (samtools view -H {input} > {output.hdo} && \
            sed '/^@SQ/{{s/SN:\([0-9XY]\)/SN:chr\\1/g;s/SN:MT/SN:chrM/g}}' {output.hdo} >  {output.hdn} && \
            samtools reheader {output.hdn} {input} > {output.bam} && \
            samtools index {output.bam})
            """
    
    rule purity_ploidy:
        input:
            normal = f"results/purity_ploidy/{config["purity_ploidy"]["normal"]}.chr.bam",
            tumor = f"results/purity_ploidy/{config["purity_ploidy"]["tumor"]}.chr.bam",
        output:
            dr = directory("results/purity_ploidy/{group}")
        conda:
            "/public/home/weiyifan/miniforge3/envs/purple"
        log:
            "logs/purity_ploidy/{group}.log"
        threads:
            5
        shell:
            """
            (export _JAVA_OPTIONS="-Xmx48g" &&\
            amber -reference ref  -reference_bam  {input.normal}  \
                -tumor tumor -tumor_bam {input.tumor} -output_dir  {output.dr}/amber \
                -loci /public/data/public_data/genomic/vcf/human/GermlineHetPon.hg38.vcf.gz \
                -ref_genome_version 38 \
                -threads {threads} && \
            cobalt -output_dir {output.dr}/cobalt -reference ref  \
                -reference_bam {input.normal} -tumor tumor -tumor_bam {input.tumor} \
                -threads {threads} \
                -gc_profile /public/home/weiyifan/software/PURPLE/copy_number/GC_profile.1000bp.37.cnp && \
            purple -output_dir {output.dr} -reference ref  -threads {threads} \
                -ensembl_data_dir /public/data/public_data/genomic/vcf/human/ensembl_data_cache \
                -amber {output.dr}/amber -cobalt {output.dr}/cobalt -tumor tumor \
                -ref_genome_version 38 -ref_genome /public/data/public_data/genomic/species/GRCh38.p14.genome.fa \
                -gc_profile /public/home/weiyifan/software/PURPLE/copy_number/GC_profile.1000bp.37.cnp > {log} 2>&1)
            """