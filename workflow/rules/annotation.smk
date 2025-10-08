rule bakta:
    """Annotate genomes with Bakta."""
    input:
        fna = "data/genomes/fna/{sample}.fna"
    output:
        gbff = "results/annotation/{sample}/{sample}.gbff",
        gff  = "results/annotation/{sample}/{sample}.gff3",
        faa  = "results/annotation/{sample}/{sample}.faa",
        tsv  = "results/annotation/{sample}/{sample}.tsv",
        json = "results/annotation/{sample}/{sample}.json",
    log:
        err = "logs/annotation/{sample}.bakta.err",
        out = "logs/annotation/{sample}.bakta.out"
    params:
        genus   = lambda wc: config["samples"][wc.sample]["genus"],
        species = lambda wc: config["samples"][wc.sample]["species"],
        strain  = lambda wc: config["samples"][wc.sample]["strain"],
        db      = config["bakta_db"]
    threads: 8
    conda: "../envs/annotation.yaml"
    shell:
        r"""
        # Create log directory if it doesn't exist
        mkdir -p logs/annotation

        # Run Bakta
        bakta --db {params.db} --threads {threads} --compliant --force \
              --genus "{params.genus}" --species "{params.species}" --strain "{params.strain}" \
              --prefix {wildcards.sample} \
              --output results/annotation/{wildcards.sample} \
              {input.fna} > {log.out} 2> {log.err}
        """