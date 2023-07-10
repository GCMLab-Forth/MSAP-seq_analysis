library(msgbsR)
library(edgeR)
library(GenomicRanges)
library(SummarizedExperiment)
library(stringr)
library(PCAtools)
library(cowplot)
library(ComplexHeatmap)

# path containing bam files to be analyzed
my_path <- "/home/msap/results/bams/" 
setwd(my_path)

se <- rawCounts(bamFilepath = my_path)

sd4 = sample.data

sd4 = data.frame(
  "Sample" = 
    c("NI_1", "NI_2", "NI_3", "NI_4",
      "INF_1", "INF_2", "INF_3", "INF_4"),
  "Treatment" =
    factor(c(0,0,0,0,1,1,1,1)))

cpmcounts <- cpm(se[rownames(ms),])
logcounts <- cpm(se,log=TRUE)

colData(se) <- DataFrame(Group =  c(rep("NI", 4), rep("INF", 4)),
                         row.names = colnames(assay(se)))

# Counts per sample plot
png("counts_per_sample_s4.png", height 
    = 4000, width = 4000, res = 480)
plotCounts(se = se, cateogory = "Group")

dev.off()

nrow(se) # 137.192 spots visible by the experiment

# Sites coordinates
cutSites <- rowRanges(se)

# Differential Methylation Analysis (DMA)
# custom threshold
ms <- diffMeth(se = se, cateogory = "Group",
               condition1 = "NI", condition2 = "INF",
               cpmThreshold = 1, thresholdSamples = 3)

a.cut = 0.05
sum(ms$FDR < a.cut) # 249 DMPs

dms = ms[ms$FDR < a.cut,]

dms.coord = as.matrix(do.call(rbind, lapply(dms[,1], function(x){unlist(str_split_fixed(x, ':|-', 4)[,-4])})))

# Lengths of contig and scaffold units
NibenLens <- read.table("genome.ctg.scf.lengths.txt",
                        header = F)

chrnames = NibenLens$V1
NibenLens = NibenLens[,2]
names(NibenLens) = chrnames

my_cuts <- GRanges(ms$site[which(ms$FDR < a.cut)])

NibenTarg = NibenLens[levels(seqnames(my_cuts))]

dmr = data.frame(seqnames(my_cuts), 
                 start(my_cuts), 
                 end(my_cuts) + 1,
                 rep("CG", times = length(my_cuts)),
                 rep(".", times = length(my_cuts)),
                 strand(my_cuts))

colnames(dmr) = c("Chr", "Start", "End", "PosType", "Score", "Strand")
write.table(dmr, file = "msgbs.p005.dms.bed", sep = '\t', quote = F, row.names = F, col.names = T)

# helper script in bash, code displayed bellow
system(command = "bash ./RMid.sh msgbs.p005.dms.bed")
# Script code
# #!/bin/bash
# 
# input=${1}
# 
# bedtools intersect -a ${1} -b my_clean_annot.gff -wao > msgbs.feature.bed
# cat msgbs.feature.bed | awk 'BEGIN{FS=OFS="\t";}  {if ($7 ~ /\./) {print $1 ":" $2 "-" $3-1 ":" $6 FS "intergenic" FS "-" FS "-" FS "." FS " "}  else {if ($15 ~ /DNA|RC|LINE|Retroposon|SINE|LTR/) {print $1 ":" $2 "-" $3-1 ":" $6 FS "repeat TE" FS $10 FS $11 FS $13 FS $15} else {print $1 ":" $2 "-" $3-1 ":" $6 FS $9 FS $10 FS $11 FS $13 FS $15}}}' > msgbs.annot.tab
# cut -f2 msgbs.annot.tab | sort | uniq -c &> fun_stats.txt
# 
# cut -f6 msgbs.annot.tab | cut -d ";" -f1 | cut -d "=" -f2 > ftype.name.tab
# cut -f1-5 msgbs.annot.tab > t && mv t msgbs.annot.tab
# paste msgbs.annot.tab ftype.name.tab > t && mv t msgbs.annot.txt
# 
# rm msgbs.feature.bed ftype.name.tab

# Microarray
(for i in `cut -f6 msgbs.annot.tab | grep "upstream" | cut -c1-22`; do grep $i ../../../../microarray.txt ;done) &> micro_matches.txt

annot.path = "msgbs.annot.txt"
msgbs.annot =  read.table(annot.path, h = F, sep = '\t')
colnames(msgbs.annot) = c("site" , "FType", "Start", "End", "Strand", "FName")

logcounts = data.frame(logcounts)
logcounts$site = rownames(logcounts)

cpmcounts = data.frame(cpmcounts)
cpmcounts$site = rownames(logcounts)

# bind with annot data 
dms_ftrs = merge(msgbs.annot, dms, by = "site")

hmd = merge(dms_ftrs, logcounts, by = "site")

print_md = merge(dms_ftrs, cpmcounts, by = "site")
print_md = print_md[,1:19]

colnames(print_md) = c(colnames(print_md)[1:11], rownames(sd4))
write.table(print_md, file = "toexcel.txt", sep = '\t', quote = F, row.names = F, col.names = T)

colnames(logcounts) = rownames(sd4)

s4 = logcounts[,1:8]
colnames(s4) = rownames(sd4)

pse_dms4 = logcounts[dms$site,][,1:8]
colnames(pse_dms4) = rownames(sd4)

PCA.CGs = 
  biplot(
    pca(s4, metadata = sd4),
    colby = "Treatment",
    title = "PCA of CGs",
    colLegendTitle = "Treatment",
    titleLabSize	= 30,
    colkey = c("0" = "light gray", "1" = "purple"), 
    pointSize = 5,
    drawConnectors = T,
    widthConnectors = 0.5,
    colConnectors = adjustcolor( "gray", alpha.f = 0.001),
    directionConnectors = "y",
    hline = 0,
    hlineCol = "dimgray",
    vline = 0,
    vlineCol = "dimgray",
    labSize = 4.5,
    min.segment.length = 1
)


PCA.DMPs = 
  biplot(pca(pse_dms4, metadata = sd4),
         colby = "Treatment",
         title = "PCA of DMPs",
         colLegendTitle = "Treatment",
         titleLabSize	= 30,
         colkey = c("0" = "light gray", "1" = "purple"), 
         pointSize = 5,
         drawConnectors = T,
         widthConnectors = 0.5,
         colConnectors = adjustcolor( "gray", alpha.f = 0.001),
         directionConnectors = "y",
         hline = 0,
         hlineCol = "dimgray",
         vline = 0,
         vlineCol = "dimgray",
         labSize = 4.5,
         min.segment.length = 1
)

png("PCAs.png", height = 4000, width = 7000, res = 480)

plot_grid(PCA.CGs, PCA.DMPs)

dev.off()


####### HEATMAP #########

s = hmd[, c(12:19)]


sLFC = hmd$logFC
ftype = hmd$FType

hm_maker = 
  function(main = "Heatmap Title", 
           fname = "heatmap.png", 
           type = F,
           rev_row = F,
           rev_col = F,
           legend=F, 
           annot_legend=F,
           sz = 1,
           name = "z-scored\nlogged cpm"
  ){
    if(type != F){
      sl = s[hmd$FType == type,]
      slLFC = sLFC[hmd$FType == type]
      ftypel = ftype[hmd$FType == type]
    }
    else{
      sl = s
      slLFC = sLFC
      ftypel = ftype
    }
    
    s.zscores = t(apply(sl, 1, function(i) ( (i-mean(as.numeric(i))) / sd(as.numeric(i)))))
    #colnames(s.zscores) = c("NI_1", "NI_2", "NI_3", "INF_1", "INF_2", "INF_3")
    colnames(s.zscores) = c("NI_1", "NI_2", "NI_3", "NI_4", "INF_1", "INF_2", "INF_3", "INF_4")
    
    q = quantile(slLFC)
    
    if (rev_col == F){
      dend_col = hclust(dist(t(s.zscores)))
    }
    else{
      dend_col = dendextend::rotate(hclust(dist(t(s.zscores))), seq(ncol(s.zscores),1))
    }
    
    if (rev_row == F){
      dend_row = hclust(dist(s.zscores))
    }
    else{
      dend_row = dendextend::rotate(hclust(dist(s.zscores)), seq(nrow(s.zscores),1))
    }
    
    png(fname, height = 5000, width = 4600, res = 450)
    #tiff(fname, height = 5000, width = 4000,res = 750)
    
    print(pheatmap(s.zscores, 
                   main = str_glue(main, " - N:", nrow(s.zscores)),
                   name = name,
                   color = brewer.pal(9, "YlOrRd"),
                   show_colnames = T,
                   show_rownames = F,
                   cluster_cols = dend_col,
                   cluster_rows = dend_row,
                   treeheight_col = 50,
                   treeheight_row = 100,
                   legend = legend,
                   fontsize = 10 * sz,
                   display_numbers= F,
                   angle_col = '0',
                   #annotation_col = data.frame("Group" = factor(rep(c("NI", "INF"), each = 4))),
                   annotation_row = data.frame("LogFC" = slLFC, "FeatureType" = ftypel),
                   annotation_colors = list(#"Group" = c("NI" = "light gray", "INF" = "purple"), 
                                            "LogFC" = diverge.color("#01665e", "#8C510A", q[1], q[5]),
                                            "FeatureType" = c("gene" = "green", 
                                                              "upstream_1000" = "blue", 
                                                              "intergenic" = "rosybrown1", 
                                                              "repeat" = "orange", 
                                                              "repeat TE" = "red")),
                   cutree_rows = 2,
                   annotation_legend = annot_legend
    ))
    dev.off()
  }


hm_maker("DMPs", "heatmap.DMPs.png", F, F, T, F, F, 2)
hm_maker("DMPs within Gene Bodies", str_glue("msgbsR.", ns, a, "DMRS.genes.png"), type = "gene", T, T)
hm_maker("DMPs within Upstream Regions", str_glue("msgbsR.", ns, a, "DMRS.upstream.png"), type = "upstream_1000", T, T)
hm_maker("DMPs within Repeats", str_glue("msgbsR.", ns, a, "DMRS.rep.png"), type = "repeat")
hm_maker("DMPs within TEs", str_glue("msgbsR.", ns, a, "DMRS.TEs.png"), type = "repeat TE", F, T)
hm_maker("DMPs within Intergenic Regions", str_glue("msgbsR.", ns, a, "DMRS.intergen.png"), type = "intergenic", T, T)