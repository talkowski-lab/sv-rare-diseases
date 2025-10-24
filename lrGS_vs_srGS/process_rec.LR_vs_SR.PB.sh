#!R
d1=read.table('LR_vs_SR.PB.bed', header=T, comment.char="", sep='\t')

anno=read.table('phase4_all_batches.annotated.bed', header=T, comment.char="", sep='\t')
colnames(anno)[4] = 'VID'
anno$VID = paste('__', anno$VID, '__', sep='')

SVID_GC = read.table('SR/PB.SVID_GC', header=T)
colnames(SVID_GC) = c('VID','GC')

SVID_filter=read.table('SR/PB.SVID_filter')
colnames(SVID_filter) = c('VID','filter')

dat=merge(d1, anno[,c('VID','PREDICTED_LOF','SVTYPE', 'samples')], by='VID')
dat=merge(dat, SVID_GC, by='VID')
dat = dat[dat$VID%in%SVID_filter[SVID_filter$filter=="PASS",]$VID,]


dat_lof = dat[!is.na(dat$PREDICTED_LOF),]
write.table(dat_lof, 'LR_vs_SR.PB.lof.bed', quote=F, sep='\t', col.names=T, row.names=F)


gene=read.table('gene_information_loeuf_pHaplo_pTriplo_pLI_4.5.txt.gz', header=T)
gene = gene[!is.na(gene$LOEUF),]
gene$LOEUF = as.double(gene$LOEUF)
constraint = gene[gene$LOEUF<.35,]$gene_name

dat_lof[,ncol(dat_lof)+1] = sapply(dat_lof$PREDICTED_LOF, function(x){min(gene[gene$gene_name%in%strsplit(as.character(x),',')[[1]], ]$LOEUF)})
dat_lof_constraint = dat_lof[dat_lof[,ncol(dat_lof)]!="Inf" & dat_lof[,ncol(dat_lof)]<.35,]
write.table(dat_lof_constraint, 'LR_vs_SR.PB.lof_constraint.bed', quote=F, sep='\t', col.names=T, row.names=F)


ad_gene = read.table('OMIM_parsed_8_April_24.AD_genes.tsv', header=T)

dat_lof[,ncol(dat_lof)+1] = sapply(dat_lof$PREDICTED_LOF, function(x){nrow(ad_gene[ad_gene$approvedGeneSymbol%in%strsplit(as.character(x),',')[[1]],]) })
dat_lof_ad = dat_lof[dat_lof[,ncol(dat_lof)]>0,]
write.table(dat_lof_ad, 'LR_vs_SR.PB.lof_AD.bed', quote=F, sep='\t', col.names=T, row.names=F)



#linux command
#to integrate lof SVs counts per genome for :

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(awk -v t="$svtype" '{if ($5==t) print}' LR_vs_SR.PB.lof.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t) print}' LR_vs_SR.PB.lof.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} all /" >> Count_per_sample.PB.lof.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(awk -v t="$svtype" '{if ($5==t && ($NF=="US" || $NF=="RM")) print}' LR_vs_SR.PB.lof.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($NF=="US" || $NF=="RM")) print}' LR_vs_SR.PB.lof.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} US_RM /" >> Count_per_sample.PB.lof.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(awk -v t="$svtype" '{if ($5==t && ($NF=="SD" || $NF=="SR")) print}' LR_vs_SR.PB.lof.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($NF=="SD" || $NF=="SR")) print}' LR_vs_SR.PB.lof.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} SD_SR /" >> Count_per_sample.PB.lof.tsv
done



#linux command
#to integrate lof_constraint SVs counts per genome for :

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(awk -v t="$svtype" '{if ($5==t) print}' LR_vs_SR.PB.lof_constraint.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t) print}' LR_vs_SR.PB.lof_constraint.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} all /" >> Count_per_sample.PB.lof_constraint.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(awk -v t="$svtype" '{if ($5==t && ($(NF-1)=="US" || $(NF-1)=="RM")) print}' LR_vs_SR.PB.lof_constraint.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($(NF-1)=="US" || $(NF-1)=="RM")) print}' LR_vs_SR.PB.lof_constraint.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} US_RM /" >> Count_per_sample.PB.lof_constraint.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(awk -v t="$svtype" '{if ($5==t && ($(NF-1)=="SD" || $(NF-1)=="SR")) print}' LR_vs_SR.PB.lof_constraint.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($(NF-1)=="SD" || $(NF-1)=="SR")) print}' LR_vs_SR.PB.lof_constraint.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} SD_SR /" >> Count_per_sample.PB.lof_constraint.tsv
done


#linux command
#to integrate lof_AD SVs counts per genome for :

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(awk -v t="$svtype" '{if ($5==t) print}' LR_vs_SR.PB.lof_AD.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t) print}' LR_vs_SR.PB.lof_AD.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} all /" >> Count_per_sample.PB.lof_AD.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(awk -v t="$svtype" '{if ($5==t && ($(NF-2)=="US" || $(NF-2)=="RM")) print}' LR_vs_SR.PB.lof_AD.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($(NF-2)=="US" || $(NF-2)=="RM")) print}' LR_vs_SR.PB.lof_AD.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} US_RM /" >> Count_per_sample.PB.lof_AD.tsv
done

for svtype in DEL DUP INS INV CPX; do
  paste \
    <(awk -v t="$svtype" '{if ($5==t && ($(NF-2)=="SD" || $(NF-2)=="SR")) print}' LR_vs_SR.PB.lof_AD.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
    <(awk -v t="$svtype" '{if ($8!="NO_OVR" && $5==t && ($(NF-2)=="SD" || $(NF-2)=="SR")) print}' LR_vs_SR.PB.lof_AD.bed \
      | awk '{print $15}' \
      | sed 's/,/\n/g' \
      | sort \
      | uniq -c \
      | awk '{sum += $1; count += 1} END {if (count > 0) print sum / count}') \
  | sed "s/^/${svtype} SD_SR /" >> Count_per_sample.PB.lof_AD.tsv
done



#!R

d1=read.table('Count_per_sample.PB.lof.tsv')
d2=read.table('Count_per_sample.PB.lof_constraint.tsv')
d3=read.table('Count_per_sample.PB.lof_AD.tsv')
dat=merge(d1, d2, by=c('V1','V2'), all=T)
dat=merge(dat, d3, by=c('V1','V2'), all=T)
dat=dat[order(dat[,1]),]
dat=dat[order(dat[,2]),]




