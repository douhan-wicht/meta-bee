# Convert xlsx to tsv
rule xlsx2tsv:
    input:
        raw = expand("data/raw/growth-data/plate-reader/{experiment}.xlsx",
                     experiment=list(config["experiments"].keys())),
        metadata = "data/raw/growth-data/plate-reader/metadata.xlsx"
    output:
        "results/preprocessing/growth-data/plate-reader/growth-data.tsv"
    log:
        out = "logs/preprocessing/xlsx2tsv.out",
        err = "logs/preprocessing/xlsx2tsv.err"
    conda:
        "../envs/preprocessing.yaml"
    params:
        raw=lambda wildcards, input: " ".join(f"--raw '{p}'" for p in input.raw)
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log.out})" "$(dirname {output})"

        python workflow/scripts/preprocessing/xlsx2tsv.py \
            {params.raw} \
            --metadata {input.metadata} \
            --output {output} \
            1> >(tee -a {log.out}) \
            2> >(tee -a {log.err} >&2)
        """

rule growth_easylinear:
    input:
        tsv = "results/preprocessing/growth-data/plate-reader/growth-data.tsv"
    output:
        summary   = "results/preprocessing/growth-data/summary/growth_easylinear.csv",
        plots_dir = directory("results/preprocessing/growth-data/summary/growth_plots")
    log:
        out = "logs/preprocessing/growth_easylinear.out",
        err = "logs/preprocessing/growth_easylinear.err"
    conda:
        "../envs/preprocessing.yaml"
    params:
        time_unit = "hours",
        h_window  = 5,
        quota     = 0.95
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log.out})" "$(dirname {output.summary})" "{output.plots_dir}"

        Rscript workflow/scripts/preprocessing/growth_easylinear.R \
          --input-tsv  {input.tsv} \
          --output-csv {output.summary} \
          --plot-dir   {output.plots_dir} \
          --time-unit  {params.time_unit} \
          --h-window   {params.h_window} \
          --quota      {params.quota} \
          1> >(tee -a {log.out}) \
          2> >(tee -a {log.err} >&2)
        """

# rule growth_amiga:
#     input:
#         tsv = "results/preprocessing/growth-data/plate-reader/growth-data.tsv"
#     output:
#         csv = "results/preprocessing/growth-data/summary/growth_amiga.csv"
#     log:
#         out = "logs/preprocessing/growth_amiga.out",
#         err = "logs/preprocessing/growth_amiga.err"
#     conda:
#         "../envs/preprocessing.yaml"
#     shell:
#         r"""
#         set -euo pipefail
#         mkdir -p "$(dirname {log.out})" "$(dirname {output.csv})"

#         python workflow/scripts/preprocessing/growth_amiga.py \
#             --input-tsv  {input.tsv} \
#             --output-csv {output.csv} \
#             1> >(tee -a {log.out}) \
#             2> >(tee -a {log.err} >&2)
#         """

rule growth_binary:
    input:
        amiga = "results/preprocessing/growth-data/summary/growth_easylinear.csv"
    output:
        "results/preprocessing/growth-data/summary/growth_binary.csv"
    log:
        out = "logs/preprocessing/growth_binary.out",
        err = "logs/preprocessing/growth_binary.err"
    conda:
        "../envs/preprocessing.yaml"
    params:
        r2_min = 0.95
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log.out})" "$(dirname {output})"

        Rscript workflow/scripts/preprocessing/growth_binary.R \
            --input  {input.amiga} \
            --output {output} \
            --r2-min {params.r2_min} \
            1> >(tee -a {log.out}) \
            2> >(tee -a {log.err} >&2)
        """

rule growth_heatmap:
    input:
        binary = "results/preprocessing/growth-data/summary/growth_binary.csv"
    output:
        png = "results/preprocessing/growth-data/summary/growth_heatmap.png",
        pdf = "results/preprocessing/growth-data/summary/growth_heatmap.pdf"
    log:
        out = "logs/preprocessing/growth_heatmap.out",
        err = "logs/preprocessing/growth_heatmap.err"
    conda:
        "../envs/preprocessing.yaml"
    params:
        # leave empty for default ordering
        species_order = "",
        media_order   = ""
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log.out})" "$(dirname {output.png})"

        Rscript workflow/scripts/preprocessing/growth_heatmap.R \
          --input       {input.binary} \
          --output-png  {output.png} \
          --output-pdf  {output.pdf} \
          {params.species_order} \
          {params.media_order} \
          1> >(tee -a {log.out}) \
          2> >(tee -a {log.err} >&2)
        """
