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

rule preprocess_models_draft:
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

rule preprocessed_extract_metrics_draft:
    """Extract metrics from PREPROCESSED draft GEMs."""
    input:
        gems = expand("results/modelProcessing/draft/{sample}/{sample}_draft_preprocessed.mat", sample=config["samples"].keys()),
        annotations = expand("results/annotation/{sample}/{sample}.txt", sample=config["samples"].keys())
    output:
        metrics = "results/modelProcessing/summary/draft_metrics.csv"
    log:
        out = "logs/modelProcessing/preprocessed_extract_metrics_draft.out",
        err = "logs/modelProcessing/preprocessed_extract_metrics_draft.err"
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

rule preprocessed_extract_metrics_relaxed:
    """Extract metrics from PREPROCESSED relaxed GEMs."""
    input:
        gems = expand("results/modelProcessing/relaxed/{sample}/{sample}_relaxed_preprocessed.mat", sample=config["samples"].keys()),
        annotations = expand("results/annotation/{sample}/{sample}.txt", sample=config["samples"].keys())
    output:
        metrics = "results/modelProcessing/summary/relaxed_metrics.csv"
    log:
        out = "logs/modelProcessing/preprocessed_extract_metrics_relaxed.out",
        err = "logs/modelProcessing/preprocessed_extract_metrics_relaxed.err"
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

rule preprocessed_extract_metrics_draft_spontaneous:
    """Extract metrics from PREPROCESSED draft_spontaneous GEMs."""
    input:
        gems = expand("results/modelProcessing/draft_spontaneous/{sample}/{sample}_draft_spontaneous_preprocessed.mat", sample=config["samples"].keys()),
        annotations = expand("results/annotation/{sample}/{sample}.txt", sample=config["samples"].keys())
    output:
        metrics = "results/modelProcessing/summary/draft_spontaneous_metrics.csv"
    log:
        out = "logs/modelProcessing/preprocessed_extract_metrics_draft_spontaneous.out",
        err = "logs/modelProcessing/preprocessed_extract_metrics_draft_spontaneous.err"
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

rule preprocessed_extract_metrics_relaxed_spontaneous:
    """Extract metrics from PREPROCESSED relaxed_spontaneous GEMs."""
    input:
        gems = expand("results/modelProcessing/relaxed_spontaneous/{sample}/{sample}_relaxed_spontaneous_preprocessed.mat", sample=config["samples"].keys()),
        annotations = expand("results/annotation/{sample}/{sample}.txt", sample=config["samples"].keys())
    output:
        metrics = "results/modelProcessing/summary/relaxed_spontaneous_metrics.csv"
    log:
        out = "logs/modelProcessing/preprocessed_extract_metrics_relaxed_spontaneous.out",
        err = "logs/modelProcessing/preprocessed_extract_metrics_relaxed_spontaneous.err"
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

rule preprocessed_model_statistics_draft:
    """Statistics for PREPROCESSED draft models."""
    input:
        metrics = "results/modelProcessing/summary/draft_metrics.csv"
    output:
        report = "results/modelProcessing/summary/draft_statistics.md"
    log:
        out = "logs/modelProcessing/preprocessed_model_statistics_draft.out",
        err = "logs/modelProcessing/preprocessed_model_statistics_draft.err"
    conda:
        "../envs/reconstruction.yaml"
    shell:
        """
        python workflow/scripts/reconstruction/model_statistics.py {input.metrics} \
            -o {output.report} \
            > {log.out} 2> {log.err}
        """

rule preprocessed_model_statistics_draft_spontaneous:
    """Statistics for PREPROCESSED draft_spontaneous models."""
    input:
        metrics = "results/modelProcessing/summary/draft_spontaneous_metrics.csv"
    output:
        report = "results/modelProcessing/summary/draft_spontaneous_statistics.md"
    log:
        out = "logs/modelProcessing/preprocessed_model_statistics_draft_spontaneous.out",
        err = "logs/modelProcessing/preprocessed_model_statistics_draft_spontaneous.err"
    conda:
        "../envs/reconstruction.yaml"
    shell:
        """
        python workflow/scripts/reconstruction/model_statistics.py {input.metrics} \
            -o {output.report} \
            > {log.out} 2> {log.err}
        """

rule preprocessed_model_statistics_relaxed:
    """Statistics for PREPROCESSED relaxed models."""
    input:
        metrics = "results/modelProcessing/summary/relaxed_metrics.csv"
    output:
        report = "results/modelProcessing/summary/relaxed_statistics.md"
    log:
        out = "logs/modelProcessing/preprocessed_model_statistics_relaxed.out",
        err = "logs/modelProcessing/preprocessed_model_statistics_relaxed.err"
    conda:
        "../envs/reconstruction.yaml"
    shell:
        """
        python workflow/scripts/reconstruction/model_statistics.py {input.metrics} \
            -o {output.report} \
            > {log.out} 2> {log.err}
        """

rule preprocessed_model_statistics_relaxed_spontaneous:
    """Statistics for PREPROCESSED relaxed_spontaneous models."""
    input:
        metrics = "results/modelProcessing/summary/relaxed_spontaneous_metrics.csv"
    output:
        report = "results/modelProcessing/summary/relaxed_spontaneous_statistics.md"
    log:
        out = "logs/modelProcessing/preprocessed_model_statistics_relaxed_spontaneous.out",
        err = "logs/modelProcessing/preprocessed_model_statistics_relaxed_spontaneous.err"
    conda:
        "../envs/reconstruction.yaml"
    shell:
        """
        python workflow/scripts/reconstruction/model_statistics.py {input.metrics} \
            -o {output.report} \
            > {log.out} 2> {log.err}
        """

rule preprocess_models_draft_spontaneous:
    """Preprocess DRAFT_SPONTANEOUS GEMs."""
    input:
        model = "results/reconstruction/draft_spontaneous/{sample}/{sample}_draft_spontaneous.mat",
        exchanges = "data/matfiles/MyExchanges.mat"
    output:
        preprocessed = "results/modelProcessing/draft_spontaneous/{sample}/{sample}_draft_spontaneous_preprocessed.mat"
    log:
        out = "logs/modelProcessing/draft_spontaneous/{sample}.out",
        err = "logs/modelProcessing/draft_spontaneous/{sample}.err"
    params:
        cobra_path  = config["cobra_path"],
        matfiles    = config["matlab_folder_path"],
        matlab_func = config["matlab_func_path"]
    threads: 4
    shell:
        r"""
        mkdir -p logs/modelProcessing/draft_spontaneous results/modelProcessing/draft_spontaneous/{wildcards.sample}

        matlab -nodisplay -nosplash -r "try; \
            addpath(genpath('{params.cobra_path}')); \
            addpath(genpath('{params.matlab_func}')); \
            preprocess_model('{params.cobra_path}', '{input.model}', '{params.matfiles}', '{output.preprocessed}'); \
        catch ME; disp(getReport(ME)); exit(1); end; exit;" \
        > {log.out} 2> {log.err}
        """

rule preprocess_models_relaxed:
    """Preprocess RELAXED GEMs."""
    input:
        model = "results/reconstruction/relaxed/{sample}/{sample}_relaxed.mat",
        exchanges = "data/matfiles/MyExchanges.mat"
    output:
        preprocessed = "results/modelProcessing/relaxed/{sample}/{sample}_relaxed_preprocessed.mat"
    log:
        out = "logs/modelProcessing/relaxed/{sample}.out",
        err = "logs/modelProcessing/relaxed/{sample}.err"
    params:
        cobra_path  = config["cobra_path"],
        matfiles    = config["matlab_folder_path"],
        matlab_func = config["matlab_func_path"]
    threads: 4
    shell:
        r"""
        mkdir -p logs/modelProcessing/relaxed results/modelProcessing/relaxed/{wildcards.sample}

        matlab -nodisplay -nosplash -r "try; \
            addpath(genpath('{params.cobra_path}')); \
            addpath(genpath('{params.matlab_func}')); \
            preprocess_model('{params.cobra_path}', '{input.model}', '{params.matfiles}', '{output.preprocessed}'); \
        catch ME; disp(getReport(ME)); exit(1); end; exit;" \
        > {log.out} 2> {log.err}
        """

rule preprocess_models_relaxed_spontaneous:
    """Preprocess RELAXED_SPONTANEOUS GEMs."""
    input:
        model = "results/reconstruction/relaxed_spontaneous/{sample}/{sample}_relaxed_spontaneous.mat",
        exchanges = "data/matfiles/MyExchanges.mat"
    output:
        preprocessed = "results/modelProcessing/relaxed_spontaneous/{sample}/{sample}_relaxed_spontaneous_preprocessed.mat"
    log:
        out = "logs/modelProcessing/relaxed_spontaneous/{sample}.out",
        err = "logs/modelProcessing/relaxed_spontaneous/{sample}.err"
    benchmark:
        "benchmarks/modelProcessing/preprocess_models_relaxed_spontaneous/{sample}.txt"
    params:
        cobra_path  = config["cobra_path"],
        matfiles    = config["matlab_folder_path"],
        matlab_func = config["matlab_func_path"]
    threads: 4
    shell:
        r"""
        mkdir -p logs/modelProcessing/relaxed_spontaneous results/modelProcessing/relaxed_spontaneous/{wildcards.sample}

        matlab -nodisplay -nosplash -r "try; \
            addpath(genpath('{params.cobra_path}')); \
            addpath(genpath('{params.matlab_func}')); \
            preprocess_model('{params.cobra_path}', '{input.model}', '{params.matfiles}', '{output.preprocessed}'); \
        catch ME; disp(getReport(ME)); exit(1); end; exit;" \
        > {log.out} 2> {log.err}
        """
