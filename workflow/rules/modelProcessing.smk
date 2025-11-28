rule add_exchange_reactions:
    """Add exchange and transport reactions for all media XLSX files."""
    input:
        matfile = "data/matfiles/GenericExchanges.mat",
        media = expand("data/media/{medium}.xlsx", medium=config["media_list"])
    output:
        updated = "data/matfiles/MyExchanges.mat"
    log:
        out = "logs/modelProcessing/add_exchange_reactions.out",
        err = "logs/modelProcessing/add_exchange_reactions.err"
    params:
        cobra_path  = config["cobra_path"],
        matlab_func = config["matlab_func_path"],
        media_args = lambda wildcards, input: ",".join(f"'{m}'" for m in input.media)
    threads: 4
    shell:
        r"""
        mkdir -p logs/modelProcessing

        matlab -nodisplay -nosplash -r "try; \
            addpath(genpath('{params.cobra_path}')); \
            addpath(genpath('{params.matlab_func}')); \
            update_exchange_file('{input.matfile}', '{output.updated}', {params.media_args}); \
        catch ME; disp(getReport(ME)); exit(1); end; exit;" \
        > {log.out} 2> {log.err}
        """

rule preprocess_models:
    """Preprocess GEMs using the custom MATLAB preprocess_model function."""
    input:
        model = "results/reconstruction/draft/{sample}/{sample}_draft.mat",
        exchanges = "data/matfiles/MyExchanges.mat"
    output:
        preprocessed = "results/modelProcessing/draft/{sample}/{sample}_draft_preprocessed.mat"
    log:
        out = "logs/modelProcessing/draft/{sample}.out",
        err = "logs/modelProcessing/draft/{sample}.err"
    params:
        cobra_path     = config["cobra_path"],
        matfiles = config["matlab_folder_path"],
        matlab_func    = config["matlab_func_path"]
    threads: 4
    shell:
        r"""
        mkdir -p logs/modelProcessing/draft results/modelProcessing/draft/{wildcards.sample}

        matlab -nodisplay -nosplash -r "try; \
            addpath(genpath('{params.cobra_path}')); \
            addpath(genpath('{params.matlab_func}')); \
            preprocess_model('{params.cobra_path}', '{input.model}', '{params.matfiles}', '{output.preprocessed}'); \
        catch ME; disp(getReport(ME)); exit(1); end; exit;" \
        > {log.out} 2> {log.err}
        """

rule extract_metrics_test_draft:
    """Extract metrics from draft GEMs."""
    input:
        gems = expand("results/modelProcessing/draft/{sample}/{sample}_draft_preprocessed.mat", sample=config["samples"].keys()),
        annotations = expand("results/annotation/{sample}/{sample}.txt", sample=config["samples"].keys())
    output:
        metrics = "results/modelProcessing/summary/draft_metrics.csv"
    log:
        out = "logs/modelProcessing/extract_metrics_draft.out",
        err = "logs/modelProcessing/extract_metrics_draft.err"
    conda:
        "../envs/reconstruction.yaml"
    shell:
        """
        python workflow/scripts/reconstruction/extract_metrics.py \
            --mat {input.gems} \
            --cds {input.annotations} \
            -o {output.metrics} \
            > {log.out} 2> {log.err}
        """

rule model_statistics_test_draft:
    """Displays reconstruction statistics (draft GEMs) and creates a markdown report to visualize them easily."""
    input:
        metrics = "results/modelProcessing/summary/draft_metrics.csv"
    output:
        report = "results/modelProcessing/summary/draft_statistics.md"
    log:
        out = "logs/modelProcessing/model_statistics_draft.out",
        err = "logs/modelProcessing/model_statistics_draft.err",
    conda:
        "../envs/reconstruction.yaml"
    shell:
        """
        python workflow/scripts/reconstruction/model_statistics.py {input.metrics} \
            -o {output.report} \
            > {log.out} 2> {log.err}
        """