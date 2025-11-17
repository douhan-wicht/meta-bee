rule draftGEMs:
    """Reconstruct draft GEMs using RAVEN's getKEGGModelForOrganism."""
    input:
        faa  = "results/annotation/{sample}/{sample}.faa"
    output:
        gem = "results/reconstruction/draft/{sample}/{sample}_draft.mat"
    log:
        out = "logs/reconstruction/draft/{sample}.out",
        err = "logs/reconstruction/draft/{sample}.err"
    params:
        raven_path = config["raven_path"],
        kegg_code = lambda wc: config["samples"][wc.sample]["kegg_code"]
    threads: 8
    shell:
        """
        # Create log and hmm directories if they don't exist
        mkdir -p data/hmm logs/reconstruction

        # Assumes RAVEN is installed and available that the path is specified in config.yaml
        matlab -nodisplay -nosplash -r "try; \
        addpath(genpath('{params.raven_path}')); \
        model = sprintf('%s_%s', '{wildcards.sample}', '{params.kegg_code}'); \
        eval([model ' = getKEGGModelForOrganism(''' '{params.kegg_code}' ''', ''' '{input.faa}' ''', ''data/hmm/prok90_kegg105'', ''results/reconstruction/draft/{wildcards.sample}'', false, false, false, false, 1e-50, 0.8, 0.3, -1);']); \
        save('{output.gem}', model); \
        catch ME; disp(getReport(ME)); exit(1); end; exit;" \
        > {log.out} 2> {log.err}
        """

rule relaxedGEMs:
    """Reconstruct relaxed GEMs using RAVEN's getKEGGModelForOrganism."""
    input:
        faa  = "results/annotation/{sample}/{sample}.faa"
    output:
        gem = "results/reconstruction/relaxed/{sample}/{sample}_relaxed.mat"
    log:
        out = "logs/reconstruction/relaxed/{sample}.out",
        err = "logs/reconstruction/relaxed/{sample}.err"
    params:
        raven_path = config["raven_path"],
        kegg_code = lambda wc: config["samples"][wc.sample]["kegg_code"]
    threads: 8
    shell:
        """
        # Create log and hmm directories if they don't exist
        mkdir -p data/hmm logs/reconstruction

        # Assumes RAVEN is installed and available that the path is specified in config.yaml
        matlab -nodisplay -nosplash -r "try; \
        addpath(genpath('{params.raven_path}')); \
        model = sprintf('%s_%s_%s', '{wildcards.sample}', '{params.kegg_code}', 'relaxed'); \
        eval([model ' = getKEGGModelForOrganism(''' '{params.kegg_code}' ''', ''' '{input.faa}' ''', ''data/hmm/prok90_kegg105'', ''results/reconstruction/relaxed/{wildcards.sample}'', false, false, false, false, 1e-10, 0.8, 0.3, -1);']); \
        save('{output.gem}', model); \
        catch ME; disp(getReport(ME)); exit(1); end; exit;" \
        > {log.out} 2> {log.err}
        """

rule extract_metrics_draft:
    """Extract metrics from draft GEMs."""
    input:
        gems = expand("results/reconstruction/draft/{sample}/{sample}_draft.mat", sample=config["samples"].keys())
    output:
        metrics = "results/reconstruction/summary/draft_metrics.csv"
    log:
        out = "logs/reconstruction/extract_metrics_draft.out",
        err = "logs/reconstruction/extract_metrics_draft.err"
    conda:
        "../envs/reconstruction.yaml"
    shell:
        """
        python workflow/scripts/reconstruction/extract_metrics.py {input.gems} \
            -o {output.metrics} \
            > {log.out} 2> {log.err}
        """

rule extract_metrics_relaxed:
    """Extract metrics from relaxed GEMs."""
    input:
        gems = expand("results/reconstruction/relaxed/{sample}/{sample}_relaxed.mat", sample=config["samples"].keys())
    output:
        metrics = "results/reconstruction/summary/relaxed_metrics.csv"
    log:
        out = "logs/reconstruction/extract_metrics_relaxed.out",
        err = "logs/reconstruction/extract_metrics_relaxed.err"
    conda:
        "../envs/reconstruction.yaml"
    shell:
        """
        python workflow/scripts/reconstruction/extract_metrics.py {input.gems} \
            -o {output.metrics} \
            > {log.out} 2> {log.err}
        """