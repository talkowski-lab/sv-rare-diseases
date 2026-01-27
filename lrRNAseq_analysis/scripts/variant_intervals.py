import polars as pl


variants = pl.read_csv('../data/lrRNA_Trios/denovo_snvs_lrRNAseq_samples.txt', separator='\t')

variants.select(
    chrom=pl.col('wgs_CHROM'),
    start=pl.col('wgs_POS') - 1,
    end=pl.col('wgs_POS') - 1 + pl.col('wgs_REF').str.len_chars(),
    id=pl.col('VarKey')
).sort(['chrom', 'start']).write_csv('../data/lrRNA_Trios/denovo_snvs.bed', separator='\t', include_header=False)

