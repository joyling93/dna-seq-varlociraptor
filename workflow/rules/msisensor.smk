if config["msisensor"]["activate"]:
    rule msisensor_scan:
        input:
            genome
        output:
            "resources/microsatellites.list"
        message:
            "Testing MSISensor scan"
        threads:
            5
        params:
            extra = ""
        log:
            "logs/msisensor_scan.log"
        conda:
            "../envs/msisensorPro.yaml"
        shell:
            """
                msisensor scan \
                -d {input} \
                -o {output} \
                > {log} 2>&1
            """
    
    rule test_msisensor_msi:
        input:
            normal = f"results/recal/{config["msisensor"]["normal"]}.bam",
            tumor = f"results/recal/{config["msisensor"]["tumor"]}.bam",
            microsat = "resources/microsatellites.list"
        output:
            "results/msisensor/{group}_dis",
            "results/msisensor/{group}_all",
            "results/msisensor/{group}_unstable"
        message:
            "Testing MSIsensor msi"
        conda:
            "../envs/msisensorPro.yaml"
        threads:
            5
        log:
            "logs/msisensor_{group}.log"
        params:
            out_prefix = "results/msisensor/{group}"
        shell:
            """
                msisensor-pro msi \
                -n {input.normal} \
                -t {input.tumor} \
                -o {params.out_prefix} \
                -d {input.microsat} \
                -b {threads} > {log} 2>&1
            """