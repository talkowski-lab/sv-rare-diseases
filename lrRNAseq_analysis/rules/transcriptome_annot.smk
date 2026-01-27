rule sqanti3_qc:
  input: 
    "data/lrRNA_Trios/transcriptomes/{sample}_R1.LRAA.gtf"
  log:
    out="logs/sqanti3_qc/{sample}.log.out",
    err="logs/sqanti3_qc/{sample}.log.err"
  params:
    output_dir=lambda w, output: os.path.dirname(str(output.iso)),
    output_prefix=lambda wildcards: wildcards.sample,
    ref_gtf=config["gencode_v48_gtf"],
    ref_fa=config["hg38_fa"]
  conda:
    "SQANTI3.env"
  output: 
    iso="results/isoform_match/sqanti3_qc/{sample}/{sample}_classification.txt",
    junc="results/isoform_match/sqanti3_qc/{sample}/{sample}_junctions.txt"
  shell: 
    """
    ~/programs/SQANTI3-5.2.2/sqanti3_qc.py {input} \
        {params.ref_gtf} \
        {params.ref_fa} \
        -d {params.output_dir} \
        -o {params.output_prefix} > {log.out} 2> {log.err}
    """

rule sqanti3_filter_ml:
  input:
    rules.sqanti3_qc.output.iso
  log:
    out="logs/sqanti3_filter_ml/{sample}.log.out",
    err="logs/sqanti3_filter_ml/{sample}.log.err"
  params:
    output_dir=lambda w, output: os.path.dirname(str(output[0])),
    output_prefix=lambda wildcards: wildcards.sample
  conda:
    "SQANTI3.env"
  output:
    "results/isoform_match/sqanti3_qc/{sample}/filter_ML/{sample}_MLresult_classification.txt"
  shell:
    """
    ~/programs/SQANTI3-5.2.2/sqanti3_filter.py ml \
        -d {params.output_dir} \
        -o {params.output_prefix} \
        {input} > {log.out} 2> {log.err}
    """

rule sqanti3_filter_rules:
  input:
    rules.sqanti3_qc.output.iso
  log:
    out="logs/sqanti3_filter_rules/{sample}.log.out",
    err="logs/sqanti3_filter_rules/{sample}.log.err"
  params:
    output_dir=lambda w, output: os.path.dirname(str(output[0])),
    output_prefix=lambda wildcards: wildcards.sample
  conda:
    "SQANTI3.env"
  output:
    "results/isoform_match/sqanti3_qc/{sample}/filter_rules/{sample}_RulesFilter_result_classification.txt"
  shell:
    """
    ~/programs/SQANTI3-5.2.2/sqanti3_filter.py rules \
        -d {params.output_dir} \
        -o {params.output_prefix} \
        {input} > {log.out} 2> {log.err}
    """

rule variant_intersect_denovo:
  input:
    junc_file=rules.sqanti3_qc.output.junc,
    var_file="data/lrRNA_Trios/denovo_snvs_lrRNAseq_samples.txt"
  params:
    prefix = lambda w, output: str(output).replace("_junc_denovo_variant_intersect.bed", "")
  conda:
    "MEI_bulk"
  output:
    "results/isoform_match/sqanti3_qc/{sample}/{sample}_junc_denovo_variant_intersect.bed"
  shell:
    """
    python scripts/denovo_variant_intersect.py \
        --sample {wildcards.sample} \
        --junc-file {input.junc_file} \
        --variant-file {input.var_file} \
        --prefix {params.prefix} 
    """

rule variant_intersect_recessive:
  input:
    junc_file=rules.sqanti3_qc.output.junc,
    var_file="data/lrRNA_Trios/rare_variants.vcf.gz"
  params:
    prefix = lambda w, output: str(output).replace("_junc_recessive_variant_intersect.bed", "")
  conda:
    "MEI_bulk"
  output:
    "results/isoform_match/sqanti3_qc/{sample}/{sample}_junc_recessive_variant_intersect.bed"
  shell:
    """
    python scripts/recessive_variant_intersect.py \
        --sample {wildcards.sample} \
        --junc-file {input.junc_file} \
        --variant-file {input.var_file} \
        --prefix {params.prefix} 
    """
