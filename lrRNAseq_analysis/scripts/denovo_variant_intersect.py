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


junc_df = pl.read_csv(args.junc_file, separator='\t')
var_df = pl.read_csv(args.variant_file, separator='\t')

_, tf_junc_name = tempfile.mkstemp(suffix=".bed")
_, tf_varfilt_name = tempfile.mkstemp(suffix=".bed")


intron_expand=20
pl.concat([
    junc_df.filter(pl.col("start_site_category") == "novel").select(
        "chrom",
        start = pl.col("genomic_start_coord") - 1 - intron_expand,
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

var_df.filter(pl.col('SAMPLE') == args.sample).select(
    chrom=pl.col('wgs_CHROM'),
    start=pl.col('wgs_POS') - 1,
    end=pl.col('wgs_POS') - 1 + pl.col('wgs_REF').str.len_chars(),
    id=pl.col('VarKey')
).sort(['chrom', 'start']).write_csv(tf_varfilt_name, separator='\t', include_header=False)

junc_bed = pbt.BedTool(tf_junc_name)
var_bed = pbt.BedTool(tf_varfilt_name)

intersected = junc_bed.intersect(var_bed, wa=True, wb=True)
intersected.saveas(f"{args.prefix}_junc_denovo_variant_intersect.bed")

