#!R
d1=read.table('SR_vs_LR.ONT.bed.gz', header=T, comment.char="", sep='\t')
d1 = d1[!is.na(d1$length) & d1$length>49,]

anno=read.table('ONT_full.annotated.bed.gz', header=T, comment.char="", sep='\t')
colnames(anno)[4] = 'VID'
anno$VID = paste('__',anno[,1],'__', anno$VID, '__', sep='')

SVID_GC = read.table('LR/ONT.SVID_GC', header=T)
colnames(SVID_GC) = c('VID','GC')
SVID_GC = unique(SVID_GC)

SVID_filter=read.table('LR/ONT.SVID_filter')
colnames(SVID_filter) = c('VID','filter')
SVID_filter = unique(SVID_filter)

dat=merge(d1, anno[,c('VID','PREDICTED_LOF','SVTYPE', 'samples')], by='VID')
dat=unique(merge(dat, unique(SVID_GC), by='VID'))
#dat=unique(merge(dat, unique(SVID_filter), by='VID'))
write.table(dat, 'SR_vs_LR.ONT.anno.bed', quote=F , sep='\t', col.names=T, row.names=F)



dat_lof = dat[!is.na(dat$PREDICTED_LOF),]
write.table(dat_lof, 'SR_vs_LR.ONT.lof.bed', quote=F, sep='\t', col.names=T, row.names=F)


gene=read.table('gene_information_loeuf_pHaplo_pTriplo_pLI_4.5.txt.gz', header=T)
gene = gene[!is.na(gene$LOEUF),]
gene$LOEUF = as.double(gene$LOEUF)
constraint = gene[gene$LOEUF<.35,]$gene_name

dat_lof[,ncol(dat_lof)+1] = sapply(dat_lof$PREDICTED_LOF, function(x){min(gene[gene$gene_name%in%strsplit(as.character(x),',')[[1]], ]$LOEUF)})
dat_lof_constraint = dat_lof[dat_lof[,ncol(dat_lof)]!="Inf" & dat_lof[,ncol(dat_lof)]<.35,]
write.table(dat_lof_constraint, 'SR_vs_LR.ONT.lof_constraint.bed', quote=F, sep='\t', col.names=T, row.names=F)


ad_gene = read.table('OMIM_parsed_8_April_24.AD_genes.tsv', header=T)

dat_lof[,ncol(dat_lof)+1] = sapply(dat_lof$PREDICTED_LOF, function(x){nrow(ad_gene[ad_gene$approvedGeneSymbol%in%strsplit(as.character(x),',')[[1]],]) })
dat_lof_ad = dat_lof[dat_lof[,ncol(dat_lof)]>0,]
write.table(dat_lof_ad, 'SR_vs_LR.ONT.lof_AD.bed', quote=F, sep='\t', col.names=T, row.names=F)



#linux command
bgzip SR_vs_LR.ONT.anno.bed
bgzip SR_vs_LR.ONT.lof_AD.bed
bgzip SR_vs_LR.ONT.lof.bed
bgzip SR_vs_LR.ONT.lof_constraint.bed



#to integrate all SVs counts per genome for :
for svtype in DEL DUP INS INV CPX; do
  paste \
    <(zcat SR_vs_LR.ONT.anno.bed.gz | awk -v t="$svtype" '{if ($5==t) print}'  \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(zcat SR_vs_LR.ONT.anno.bed.gz | awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t) print}'  \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} all /" >> Count_per_sample.LR_ONT.all.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(zcat SR_vs_LR.ONT.anno.bed.gz | awk -v t="$svtype" '{if ($5==t && ($NF=="US" || $NF=="RM")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(zcat SR_vs_LR.ONT.anno.bed.gz | awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($NF=="US" || $NF=="RM")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} US_RM /" >> Count_per_sample.LR_ONT.all.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(zcat SR_vs_LR.ONT.anno.bed.gz | awk -v t="$svtype" '{if ($5==t && ($NF=="SD" || $NF=="SR")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(zcat SR_vs_LR.ONT.anno.bed.gz | awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($NF=="SD" || $NF=="SR")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} SD_SR /" >> Count_per_sample.LR_ONT.all.tsv
done


#to integrate lof SVs counts per genome for :
for svtype in DEL DUP INS INV CPX; do
  paste \
    <(zcat SR_vs_LR.ONT.lof.bed.gz | awk -v t="$svtype" '{if ($5==t) print}'  \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(zcat SR_vs_LR.ONT.lof.bed.gz | awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t) print}'  \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} all /" >> Count_per_sample.LR_ONT.lof.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(zcat SR_vs_LR.ONT.lof.bed.gz | awk -v t="$svtype" '{if ($5==t && ($NF=="US" || $NF=="RM")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(zcat SR_vs_LR.ONT.lof.bed.gz | awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($NF=="US" || $NF=="RM")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} US_RM /" >> Count_per_sample.LR_ONT.lof.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(zcat SR_vs_LR.ONT.lof.bed.gz | awk -v t="$svtype" '{if ($5==t && ($NF=="SD" || $NF=="SR")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(zcat SR_vs_LR.ONT.lof.bed.gz | awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($NF=="SD" || $NF=="SR")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} SD_SR /" >> Count_per_sample.LR_ONT.lof.tsv
done


#linux command
#to integrate lof_constraint SVs counts per genome for :
for svtype in DEL DUP INS INV CPX; do
  paste \
    <(zcat SR_vs_LR.ONT.lof_constraint.bed.gz | awk -v t="$svtype" '{if ($5==t) print}'  \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(zcat SR_vs_LR.ONT.lof_constraint.bed.gz | awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t) print}'  \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} all /" >> Count_per_sample.LR_ONT.lof_constraint.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(zcat SR_vs_LR.ONT.lof_constraint.bed.gz | awk -v t="$svtype" '{if ($5==t && ($(NF-1)=="US" || $(NF-1)=="RM")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(zcat SR_vs_LR.ONT.lof_constraint.bed.gz | awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($(NF-1)=="US" || $(NF-1)=="RM")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} US_RM /" >> Count_per_sample.LR_ONT.lof_constraint.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(zcat SR_vs_LR.ONT.lof_constraint.bed.gz | awk -v t="$svtype" '{if ($5==t && ($(NF-1)=="SD" || $(NF-1)=="SR")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(zcat SR_vs_LR.ONT.lof_constraint.bed.gz | awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($(NF-1)=="SD" || $(NF-1)=="SR")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} SD_SR /" >> Count_per_sample.LR_ONT.lof_constraint.tsv
done


#linux command
#to integrate lof_AD SVs counts per genome for :
for svtype in DEL DUP INS INV CPX; do
  paste \
    <(zcat SR_vs_LR.ONT.lof_AD.bed.gz | awk -v t="$svtype" '{if ($5==t) print}'  \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(zcat SR_vs_LR.ONT.lof_AD.bed.gz | awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t) print}'  \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} all /" >> Count_per_sample.LR_ONT.lof_AD.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(zcat SR_vs_LR.ONT.lof_AD.bed.gz | awk -v t="$svtype" '{if ($5==t && ($(NF-2)=="US" || $(NF-2)=="RM")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(zcat SR_vs_LR.ONT.lof_AD.bed.gz | awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($(NF-2)=="US" || $(NF-2)=="RM")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} US_RM /" >> Count_per_sample.LR_ONT.lof_AD.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(zcat SR_vs_LR.ONT.lof_AD.bed.gz | awk -v t="$svtype" '{if ($5==t && ($(NF-2)=="SD" || $(NF-2)=="SR")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(zcat SR_vs_LR.ONT.lof_AD.bed.gz | awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($(NF-2)=="SD" || $(NF-2)=="SR")) print}' \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} SD_SR /" >> Count_per_sample.LR_ONT.lof_AD.tsv
done





#!R

readin_lr_sr_comp_stat<-function(file_name){
  d1=read.table(file_name, sep='\t')
  d1[,ncol(d1) +1] = sapply(d1[,1], function(x){strsplit(as.character(x),' ')[[1]][1]})
  d1[,ncol(d1) +1] = sapply(d1[,1], function(x){strsplit(as.character(x),' ')[[1]][2]})
  d1[,ncol(d1) +1] = sapply(d1[,1], function(x){strsplit(as.character(x),' ')[[1]][3]})
  out = d1[,c(3:5,2)]
  out[is.na(out)] = 0
  colnames(out)=c('svtype','gc','all','ovr')

  out$ovr = as.double(out$ovr)
  out$all = as.double(out$all)
  out[,ncol(out)+1] = out$ovr/out$all
  colnames(out)[ncol(out)] = 'ovr_rate'
  
  out$ovr = as.integer(out$ovr)
  out$all = as.integer(out$all)

  return(out)
}

d0 = readin_lr_sr_comp_stat('Count_per_sample.LR_ONT.all.tsv')
d1 = readin_lr_sr_comp_stat('Count_per_sample.LR_ONT.lof.tsv')
d2 = readin_lr_sr_comp_stat('Count_per_sample.LR_ONT.lof_constraint.tsv')
d3 = readin_lr_sr_comp_stat('Count_per_sample.LR_ONT.lof_AD.tsv')

dat=merge(d0, d1, by=c('svtype','gc'), all=T)
dat=merge(dat, d2, by=c('svtype','gc'), all=T)
dat=merge(dat, d3, by=c('svtype','gc'), all=T)
dat=dat[order(dat[,1]),]
dat=dat[order(dat[,2]),]
dat[is.na(dat)] = 0
write.table(dat, 'Count_per_sample.LR_ONT.integrated.tsv', quote=F, sep='\t', col.names=T, row.names=F)












