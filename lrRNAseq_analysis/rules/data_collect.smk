rule download_transcriptomes:
  params:
    folder="gs://fc-f3a757f6-a0e3-4475-b972-39bb3ff1f948/downstreamProcessingOuts/LRAA_refGuidedID"
  output:
    gtf="data/lrRNA_Trios/transcriptomes/{sample}_R1.LRAA.gtf",
    quant="data/lrRNA_Trios/transcriptomes/{sample}_R1.LRAA.quant.expr"
  shell:
    """
    mkdir -p data/lrRNA_Trios/transcriptomes
    gcloud storage cp {params.folder}/{wildcards.sample}_R1.LRAA.gtf {output.gtf}
    sed -i '/^$/d' {output.gtf}  # Remove empty lines
    gcloud storage cp {params.folder}/{wildcards.sample}_R1.LRAA.quant.expr {output.quant}
    """

