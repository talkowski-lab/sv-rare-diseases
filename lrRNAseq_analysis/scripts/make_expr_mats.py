from pathlib import Path
import polars as pl
import polars.selectors as csc
import subprocess


samples = sorted([
    f.name for f in Path("../results/isoform_match/sqanti3_qc").glob("RGP*")
])

transcript_counts = []
gene_counts = []


for s in samples:
    matchup_df = pl.read_csv(f"../results/isoform_match/sqanti3_qc/{s}/{s}_classification.txt", separator="\t", null_values="NA")

    quant_df = pl.read_csv(f"../data/lrRNA_Trios/transcriptomes/{s}_R1.LRAA.quant.expr", separator="\t")

    transcript_count_df = matchup_df[['isoform', 'associated_gene', 'associated_transcript']].join(
        quant_df[['transcript_id', 'all_reads', 'TPM']],
        left_on="isoform",
        right_on="transcript_id"
    ).filter(
        ~pl.col("associated_gene").str.starts_with("novel"),
        ~pl.col("associated_transcript").str.starts_with("novel"),
        ~pl.col("associated_gene").str.contains("_")
    )
    transcript_counts.append(
        transcript_count_df.group_by(["associated_gene", "associated_transcript"]).agg(
            pl.col("all_reads").sum().alias(f"{s}_count"),
            pl.col("uniq_reads").sum().alias(f"{s}_uniq_count"),
            pl.col("TPM").sum().alias(f"{s}_TPM")
        ).rename({"associated_gene": "gene_id", "associated_transcript": "transcript_id"})
    )

    gene_counts.append(
        transcript_count_df.group_by("associated_gene").agg(
            pl.col("all_reads").sum().alias(f"{s}_count"),
            pl.col("uniq_reads").sum().alias(f"{s}_uniq_count"),
            pl.col("TPM").sum().alias(f"{s}_TPM")
        ).rename({"associated_gene": "gene_id"})
    )

# Combined into single dataframe

def combined_dfs(dfs, id_cols, null_fill=0):
    if type(id_cols) is not list:
        id_cols = [id_cols]
    all_ids = pl.concat([df[id_cols] for df in dfs], how="vertical").unique().sort(by=id_cols)
    return pl.concat(
        [all_ids] +
        [all_ids.join(df, on=id_cols, how="left", coalesce=True).select(pl.exclude(id_cols)) for df in dfs]
    , how='horizontal').fill_null(null_fill)

all_genes = combined_dfs(gene_counts, id_cols = "gene_id")
all_transcripts = combined_dfs(transcript_counts, id_cols=["gene_id", "transcript_id"])

gene_counts = all_genes.select('gene_id', csc.ends_with("_count"))
gene_counts.columns = [c.replace("_count", "") for c in gene_counts.columns]

gene_uniq_counts = all_genes.select('gene_id', csc.ends_with("_uniq_count"))
gene_uniq_counts.columns = [c.replace("_uniq_count", "") for c in gene_uniq_counts.columns]

gene_TPM = all_genes.select('gene_id', csc.ends_with("_TPM"))
gene_TPM.columns = [c.replace("_TPM", "") for c in gene_TPM.columns]

transcript_counts = all_transcripts.select("gene_id", "transcript_id", csc.ends_with("_count"))
transcript_counts.columns = [c.replace("_count", "") for c in transcript_counts.columns]

transcript_TPM = all_transcripts.select('gene_id', 'transcript_id', csc.ends_with("_TPM"))
transcript_TPM.columns = [c.replace("_TPM", "") for c in transcript_TPM.columns]

transcript_uniq_counts = all_transcripts.select("gene_id", "transcript_id", csc.ends_with("_uniq_count"))
transcript_uniq_counts.columns = [c.replace("_uniq_count", "") for c in transcript_uniq_counts.columns]


gene_counts.write_csv("../results/isoform_match/gene_counts_unfiltered.txt", separator="\t")
gene_uniq_counts.write_csv("../results/isoform_match/gene_uniq_counts_unfiltered.txt", separator="\t")
gene_TPM.write_csv("../results/isoform_match/gene_TPM_unfiltered.txt", separator="\t")
 
transcript_counts.write_csv("../results/isoform_match/transcript_counts_unfiltered.txt", separator="\t")
transcript_uniq_counts.write_csv("../results/isoform_match/transcript_uniq_counts_unfiltered.txt", separator="\t")
transcript_TPM.write_csv("../results/isoform_match/transcript_TPM_unfiltered.txt", separator="\t")

for f in Path("../results/isoform_match").glob("*_unfiltered.txt"):
    subprocess.run(["gzip", str(f)])

