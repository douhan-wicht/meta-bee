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

rule growthrates_easylinear:
    input:
        tsv = "results/preprocessing/growth-data/plate-reader/growth-data.tsv"
    output:
        summary   = "results/preprocessing/growth-data/summary/growth_rates_easylinear.csv",
        plots_dir = directory("results/preprocessing/growth-data/summary/growth_plots")
    log:
        out = "logs/preprocessing/growthrates_easylinear.out",
        err = "logs/preprocessing/growthrates_easylinear.err"
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

        Rscript workflow/scripts/preprocessing/growthrates_easylinear.R \
          --input-tsv  {input.tsv} \
          --output-csv {output.summary} \
          --plot-dir   {output.plots_dir} \
          --time-unit  {params.time_unit} \
          --h-window   {params.h_window} \
          --quota      {params.quota} \
          1> >(tee -a {log.out}) \
          2> >(tee -a {log.err} >&2)
        """
