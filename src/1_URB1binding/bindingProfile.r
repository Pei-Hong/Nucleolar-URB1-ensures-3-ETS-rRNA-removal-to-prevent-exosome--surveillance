suppressMessages({library(pheatmap)
                library(getopt)})

args <- commandArgs(trailingOnly = TRUE)
bgfile = args[1]
file = read.table(bgfile, header=F, stringsAsFactors=F)

file.mat = as.matrix(file[, c("V4", "V6")])
rownames(file.mat) = file$V2
colnames(file.mat) = c("21471", "21472")

pdf("combine2_region.pdf", width= 12, height = 2)
p<- pheatmap(t(file.mat[12850:12994,]), 
    show_rownames = T, show_colnames = F, cluster_cols = F, cluster_rows = F, 
    color =  colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(100),
     legend_breaks=seq(0,0.022,0.005), legend_labels = seq(0,0.022,0.005), breaks = seq(0,0.022,0.00022) )
dev.off()