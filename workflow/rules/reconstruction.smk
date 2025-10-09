rule getKEGGModelForOrganism:
    """Reconstruct draft GEMs using RAVEN's getKEGGModelForOrganism."""
    input:
        faa  = "results/annotation/{sample}/{sample}.faa"
    output:
        gem = "results/reconstruction/{sample}/{sample}_draftGEM.mat"
    log:
        out = "logs/reconstruction/{sample}.out",
        err = "logs/reconstruction/{sample}.err"
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
        eval([model ' = getKEGGModelForOrganism(''' '{params.kegg_code}' ''', ''' '{input.faa}' ''', ''data/hmm/prok90_kegg105'', ''results/reconstruction/{wildcards.sample}'', false, false, false, false, 1e-50, 0.8, 0.3, -1);']); \
        save('{output.gem}', model); \
        catch ME; disp(getReport(ME)); exit(1); end; exit;" \
        > {log.out} 2> {log.err}
        """