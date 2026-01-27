import tempfile
import argparse
import polars as pl
import pybedtools as pbt


parser = argparse.ArgumentParser()
parser.add_argument('--sample')
parser.add_argument('--junc-file')
parser.add_argument('--variant-file')
parser.add_argument('--prefix')
args = parser.parse_args()


# args = argparse.Namespace()
# args.sample = "RGP_657_3"
# args.junc_file = "../results/isoform_match/sqanti3_qc/RGP_657_3/RGP_657_3_junctions.txt"
# args.variant_file = "../data/lrRNA_Trios/rare_variants.vcf"
# args.prefix = "../results/isoform_match/sqanti3_qc/RGP_657_3/RGP_657_3"

junc_df = pl.read_csv(args.junc_file, separator='\t')

tf_junc_name = tempfile.mktemp(suffix=".bed")



intron_expand=20
pl.concat([
    junc_df.filter(pl.col("start_site_category") == "novel").select(
        "chrom",
        start = pl.col("genomic_start_coord") - 1 + intron_expand,
        end=pl.col("genomic_start_coord") + 1 + intron_expand,
        id=pl.concat_str(["isoform", pl.col("junction_number").str.replace("junction_", ""), pl.lit("start")], separator=":")
    ),
    junc_df.filter(pl.col("end_site_category") == "novel").select(
        "chrom",
        start = pl.col("genomic_end_coord") - 2 - intron_expand,
        end=pl.col("genomic_end_coord") + intron_expand,
        id=pl.concat_str(["isoform", pl.col("junction_number").str.replace("junction_", ""), pl.lit("end")], separator=":")
    ),

]).sort("chrom", "start").write_csv(tf_junc_name, separator='\t', include_header=False)

junc_bed = pbt.BedTool(tf_junc_name)
var_bed = pbt.BedTool(args.variant_file)

intersected = junc_bed.intersect(var_bed, wa=True, wb=True)

intersected_df = pl.read_csv(intersected.fn, separator='\t', has_header=False).drop("column_5", "column_10", "column_11", "column_13")
intersected_df.columns = ["chrom", "start", "end", "junction_id", "POS", "ID", "REF", "ALT", "INFO"]
intersected_df.filter(pl.col("ID").str.contains(args.sample)).with_columns(
    var_start=pl.col("POS") - 1,
    var_end=pl.when(pl.col("INFO").str.contains("END=")).then(
        pl.col("INFO").str.extract(r"END=(\d+)", 1).cast(pl.Int64)
    ).otherwise(pl.col("POS") + pl.col("REF").str.len_chars() - 1)
).drop("POS", "REF", "ALT", "INFO").write_csv(f"{args.prefix}_junc_recessive_variant_intersect.bed", separator="\t", include_header=False)

