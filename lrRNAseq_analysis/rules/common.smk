with open('data/lrRNA_Trios/families.txt', 'r') as inf:
  families = [line.strip() for line in inf]

ch_samples = [f"{fam}_3" for fam in families]
samples = [f"{fam}_{i}" for i in [1,2,3] for fam in families] 

def get_final_outputs():
  o = []
  o.extend(expand(rules.download_transcriptomes.output.gtf, sample=samples))
  o.extend(expand(rules.sqanti3_filter_ml.output[0], sample=samples))
  o.extend(expand(rules.sqanti3_filter_rules.output[0], sample=samples))
  o.extend(expand(rules.variant_intersect_denovo.output, sample=ch_samples))
  o.extend(expand(rules.variant_intersect_recessive.output, sample=ch_samples))
  return o
