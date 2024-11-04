if config["driver_gene"]["activate"]:
    rule driver_gene:
        input:
            "results/tables/{group}.{event}.{calling_type}.fdr-controlled.tsv",
        output:
            "results/driver_gene/{group}.{event}.{calling_type}.driver_gene.xls",
        log:
            "logs/driver_gene/{group}.{event}.{calling_type}.log"
        script:
            "../scripts/driver.py"