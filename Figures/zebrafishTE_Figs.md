### R code for figure generation ###

Bo Li, Puwen Tan and Yunhao Wang



```R
## Fig.1g

library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(gplots)
library (vegan)
library(viridis)
#library(ComplexHeatmap)

paletteLength <- 50
myColor2 <- colorRampPalette(c("#EBEBEB","white", "#FFE949"))(paletteLength)
pdf("ERV_core_domain.sum.Heatmap.2.pdf",w=6,h=6)
#annotation_row<-read.table("DNA_subfamily_annotation",header=T,sep="\t",,row.names=1)
#annotation_row <- as.data.frame(annotation_row)
mydata <- read.table("ERV_core_domain.sum.4R", header=T, sep="\t")
pheatmap(mydata, scale = "none", color = myColor2, border_color = "#EBEBEB", cluster_cols=FALSE,cluster_rows=FALSE)
dev.off()

paletteLength <- 50
myColor2 <- colorRampPalette(c("#EBEBEB","white", "#FFCB49"))(paletteLength)
pdf("LTR_core_domain.sum.4R.Heatmap.2.pdf",w=6,h=6)
#annotation_row<-read.table("DNA_subfamily_annotation",header=T,sep="\t",,row.names=1)
#annotation_row <- as.data.frame(annotation_row)
mydata <- read.table("LTR_core_domain.sum.4R", header=T, sep="\t")
pheatmap(mydata, scale = "none", color = myColor2, border_color = "#EBEBEB",cluster_cols=FALSE,cluster_rows=FALSE)
dev.off()

paletteLength <- 50
myColor2 <- colorRampPalette(c("#EBEBEB","white", "#3D8AD5"))(paletteLength)
pdf("DNA_core_domain.sum.4R.Heatmap.2.pdf",w=6,h=6)
#annotation_row<-read.table("DNA_subfamily_annotation",header=T,sep="\t",,row.names=1)
#annotation_row <- as.data.frame(annotation_row)
mydata <- read.table("DNA_core_domain.sum.4R", header=T, sep="\t")
pheatmap(mydata, scale = "none", color = myColor2, border_color = "#EBEBEB",cluster_cols=FALSE,cluster_rows=FALSE)
dev.off()


paletteLength <- 50
myColor2 <- colorRampPalette(c("#EBEBEB","white", "#F18F48"))(paletteLength)
pdf("LINE_core_domain.sum.4R.Heatmap.2.pdf",w=6,h=6)
#annotation_row<-read.table("DNA_subfamily_annotation",header=T,sep="\t",,row.names=1)
#annotation_row <- as.data.frame(annotation_row)
mydata <- read.table("LINE_core_domain.sum.4R", header=T, sep="\t")
pheatmap(mydata, scale = "none", color = myColor2, border_color = "#EBEBEB", cluster_cols=FALSE,cluster_rows=FALSE)
dev.off()
```


```R
## Fig2b
#install.packages("ggcorrplot")
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggcorrplot)

data <- read.table("TE_trans.final.cpm.noRep.matrix", header=T, row.names=1, com='', quote='', check.names=F, sep="\t")
cormat <- round(cor(data, method="pearson"),2)
melted_cormat <- melt(cormat)

pdf("TE_trans.final.cpm.noRep.matrix.pearson.lower.1.pdf", pointsize=5)

ggcorrplot(cormat, hc.order = FALSE, type = "lower",
            outline.col = "white", lab=TRUE, lab_col="black")+
  #scale_fill_gradientn(colours = colorspace::diverge_hcl(7))+
  #scale_fill_gradient(low = "white", high = "#A9D18E")+
  #scale_fill_gradient2(low = "grey",mid = "white", high ="#A9D18E",midpoint=0.2 )+
  scale_fill_gradientn(colours = colorspace::sequential_hcl(7, h = c(260, 60), c = 60, l = c(40, 95), power = 1))+
  #scale_fill_gradient2(low = "#A9D18E",mid = "white", high = "grey",midpoint=0.2 )+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 14, hjust = 1))+theme(axis.text.y = element_text(size = 14, hjust = 1))
                                                      
dev.off()


data <- read.table("Gene_trans.final.cpm.noRep", header=T, row.names=1, com='', quote='', check.names=F, sep="\t")
cormat <- round(cor(data, method="pearson"),2)
melted_cormat <- melt(cormat)

pdf("Gene_trans.final.cpm.noRep.matrix.pearson.lower.1.pdf", pointsize=5)

ggcorrplot(cormat, hc.order = FALSE, type = "lower",
            outline.col = "white", lab=TRUE, lab_col="black")+
  #scale_fill_gradientn(colours = colorspace::diverge_hcl(7))+
  #scale_fill_gradient2(low = "grey",mid = "white", high ="#A9D18E",midpoint=0.05 )+
  scale_fill_gradientn(colours = colorspace::sequential_hcl(7, h = c(260, 60), c = 60, l = c(40, 95), power = 1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 14, hjust = 1))+theme(axis.text.y = element_text(size = 14, hjust = 1))
                                                      
dev.off()

data <- read.table("TE-Gene_trans.final.cpm.noRep", header=T, row.names=1, com='', quote='', check.names=F, sep="\t")
cormat <- round(cor(data, method="pearson"),2)
melted_cormat <- melt(cormat)

pdf("TE-Gene_trans.final.cpm.noRep.matrix.pearson.lower.1.pdf", pointsize=5)

ggcorrplot(cormat, hc.order = FALSE, type = "lower",
            outline.col = "white", lab=TRUE, lab_col="black")+
  scale_fill_gradientn(colours = colorspace::sequential_hcl(7, h = c(260, 60), c = 60, l = c(40, 95), power = 1))+
  #scale_fill_gradient2(low = "grey",high ="#0DD10E",midpoint=0.05 )+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 14, hjust = 1))+theme(axis.text.y = element_text(size = 14, hjust = 1))
                                                      
dev.off()


```


```R
## Fig.2d-f and Supplementary Fig.8
library(ggplot2)
#library(tidyverse)


# Set the output PDF file
pdf("Aman_4hpf_vs_WT_4hpf.scatter.pdf", w = 8, h = 6)

# Read your data from the file
data <- read.table("Aman_4hpf_vs_WT_4hpf.scatter.txt", header = TRUE, sep = "\t")

data1 <- data[data$WT >= 1 & data$Class != "Maternal_gene", ]
data2 <-rbind (data[data$Class == "Maternal_gene", ], data1)

# Create the scatter plot with different alpha values
ggplot(data2, aes(x = WT, y = Mutant, color = Class)) +
  geom_point(aes(alpha = Class), size = 4) +
  scale_color_manual(values = c('grey','#B02418','#000000')) +
  scale_alpha_manual(values = c(0.1, 0.8, 0.8)) +
  scale_x_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  scale_y_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  geom_abline(lty = 1) +
  #scale_alpha_manual(values = c('grey' = 0.1, '#E69F00' = 1, '#000000' = 1)) +  # Map alpha values to Class
  theme_classic()+ coord_fixed()

# Save the plot to the PDF file
dev.off()


# Set the output PDF file
pdf("Aman_6hpf_vs_WT_6hpf.scatter.pdf", w = 8, h = 6)

# Read your data from the file
data <- read.table("Aman_6hpf_vs_WT_6hpf.scatter.txt", header = TRUE, sep = "\t")

# Create the scatter plot with different alpha values
ggplot(data, aes(x = WT, y = Mutant, color = Class)) +
  geom_point(aes(alpha = Class), size = 4) +
  scale_color_manual(values = c('grey','#B02418','#000000')) +
  scale_alpha_manual(values = c(0.1, 0.8, 0.8)) +
  scale_x_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  scale_y_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  geom_abline(lty = 1) +
  #scale_alpha_manual(values = c('grey' = 0.1, '#E69F00' = 1, '#000000' = 1)) +  # Map alpha values to Class
  theme_classic()+ coord_fixed()

# Save the plot to the PDF file
dev.off()


# Set the output PDF file
pdf("CHX_4hpf_vs_WT_4hpf.scatter.pdf", w = 8, h = 6)

# Read your data from the file
data <- read.table("CHX_4hpf_vs_WT_4hpf.scatter.txt", header = TRUE, sep = "\t")

data1 <- data[data$WT >= 1 & data$Class != "Maternal_gene", ]
data2 <-rbind (data[data$Class == "Maternal_gene", ], data1)

# Create the scatter plot with different alpha values
ggplot(data2, aes(x = WT, y = Mutant, color = Class)) +
  geom_point(aes(alpha = Class), size = 4) +
  scale_color_manual(values = c('grey','#B02418','#000000')) +
  scale_alpha_manual(values = c(0.1, 0.8, 1)) +
  scale_x_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  scale_y_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  geom_abline(lty = 1) +
  #scale_alpha_manual(values = c('grey' = 0.1, '#E69F00' = 1, '#000000' = 1)) +  # Map alpha values to Class
  theme_classic()+ coord_fixed()

# Save the plot to the PDF file
dev.off()


# Set the output PDF file
pdf("CHX_6hpf_vs_WT_6hpf.scatter.pdf", w = 8, h = 6)

# Read your data from the file
data <- read.table("CHX_6hpf_vs_WT_6hpf.scatter.txt", header = TRUE, sep = "\t")

data1 <- data[data$WT >= 1 & data$Class != "Maternal_gene", ]
data2 <-rbind (data[data$Class == "Maternal_gene", ], data1)

# Create the scatter plot with different alpha values
ggplot(data2, aes(x = WT, y = Mutant, color = Class)) +
  geom_point(aes(alpha = Class), size = 4) +
  scale_color_manual(values = c('grey','#B02418','#000000')) +
  scale_alpha_manual(values = c(0.1, 0.8, 1)) +
  scale_x_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  scale_y_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  geom_abline(lty = 1) +
  #scale_alpha_manual(values = c('grey' = 0.1, '#E69F00' = 1, '#000000' = 1)) +  # Map alpha values to Class
  theme_classic()+ coord_fixed()

# Save the plot to the PDF file
dev.off()


# Set the output PDF file
pdf("tri_res_4hpf_vs_WT_4hpf.scatter.pdf", w = 8, h = 6)

# Read your data from the file
data <- read.table("tri_res_4hpf_vs_WT_4hpf.scatter.txt", header = TRUE, sep = "\t")

# Create the scatter plot with different alpha values
ggplot(data, aes(x = WT, y = Mutant, color = Class)) +
  geom_point(aes(alpha = Class), size = 4) +
  scale_color_manual(values = c('grey','#B02418','#000000')) +
  scale_alpha_manual(values = c(0.1, 0.8, 1)) +
  scale_x_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  scale_y_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  geom_abline(lty = 1) +
  #scale_alpha_manual(values = c('grey' = 0.1, '#E69F00' = 1, '#000000' = 1)) +  # Map alpha values to Class
  theme_classic()+ coord_fixed()

# Save the plot to the PDF file
dev.off()

# Set the output PDF file
pdf("tri_res_6hpf_vs_WT_6hpf.scatter.pdf", w = 8, h = 6)

# Read your data from the file
data <- read.table("tri_res_6hpf_vs_WT_6hpf.scatter.txt", header = TRUE, sep = "\t")

# Create the scatter plot with different alpha values
ggplot(data, aes(x = WT, y = Mutant, color = Class)) +
  geom_point(aes(alpha = Class), size = 4) +
  scale_color_manual(values = c('grey','#B02418','#000000')) +
  scale_alpha_manual(values = c(0.1, 0.8, 1)) +
  scale_x_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  scale_y_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  geom_abline(lty = 1) +
  #scale_alpha_manual(values = c('grey' = 0.1, '#E69F00' = 1, '#000000' = 1)) +  # Map alpha values to Class
  theme_classic()+ coord_fixed()

# Save the plot to the PDF file
dev.off()

# Set the output PDF file
pdf("tri_res_8hpf_vs_WT_8hpf.scatter.pdf", w = 8, h = 6)

# Read your data from the file
data <- read.table("tri_res_8hpf_vs_WT_8hpf.scatter.txt", header = TRUE, sep = "\t")

# Create the scatter plot with different alpha values
ggplot(data, aes(x = WT, y = Mutant, color = Class)) +
  geom_point(aes(alpha = Class), size = 4) +
  scale_color_manual(values = c('grey','#B02418','#000000')) +
  scale_alpha_manual(values = c(0.1, 0.8, 1)) +
  scale_x_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  scale_y_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  geom_abline(lty = 1) +
  #scale_alpha_manual(values = c('grey' = 0.1, '#E69F00' = 1, '#000000' = 1)) +  # Map alpha values to Class
  theme_classic()+ coord_fixed()

# Save the plot to the PDF file
dev.off()


# Set the output PDF file
pdf("tri_4hpf_vs_WT_4hpf.scatter.pdf", w = 8, h = 6)

# Read your data from the file
data <- read.table("tri_4hpf_vs_WT_4hpf.scatter.txt", header = TRUE, sep = "\t")

# Create the scatter plot with different alpha values
ggplot(data, aes(x = WT, y = Mutant, color = Class)) +
  geom_point(aes(alpha = Class), size = 4) +
  scale_color_manual(values = c('grey','#B02418','#000000')) +
  scale_alpha_manual(values = c(0.1, 0.8, 1)) +
  scale_x_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  scale_y_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  geom_abline(lty = 1) +
  #scale_alpha_manual(values = c('grey' = 0.1, '#E69F00' = 1, '#000000' = 1)) +  # Map alpha values to Class
  theme_classic()+ coord_fixed()

# Save the plot to the PDF file
dev.off()


# Set the output PDF file
pdf("tri_6hpf_vs_WT_6hpf.scatter.pdf", w = 8, h = 6)

# Read your data from the file
data <- read.table("tri_6hpf_vs_WT_6hpf.scatter.txt", header = TRUE, sep = "\t")

# Create the scatter plot with different alpha values
ggplot(data, aes(x = WT, y = Mutant, color = Class)) +
  geom_point(aes(alpha = Class), size = 4) +
  scale_color_manual(values = c('grey','#B02418','#000000')) +
  scale_alpha_manual(values = c(0.1, 0.8, 1)) +
  scale_x_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  scale_y_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  geom_abline(lty = 1) +
  #scale_alpha_manual(values = c('grey' = 0.1, '#E69F00' = 1, '#000000' = 1)) +  # Map alpha values to Class
  theme_classic()+ coord_fixed()

# Save the plot to the PDF file
dev.off()

# Set the output PDF file
pdf("tri_8hpf_vs_WT_8hpf.scatter.pdf", w = 8, h = 6)

# Read your data from the file
data <- read.table("tri_8hpf_vs_WT_8hpf.scatter.txt", header = TRUE, sep = "\t")

# Create the scatter plot with different alpha values
ggplot(data, aes(x = WT, y = Mutant, color = Class)) +
  geom_point(aes(alpha = Class), size = 4) +
  scale_color_manual(values = c('grey','#B02418','#000000')) +
  scale_alpha_manual(values = c(0.1, 0.8, 1)) +
  scale_x_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  scale_y_continuous( breaks = c(0, 2,4,6,8,10,12,14))+
  geom_abline(lty = 1) +
  #scale_alpha_manual(values = c('grey' = 0.1, '#E69F00' = 1, '#000000' = 1)) +  # Map alpha values to Class
  theme_classic()+ coord_fixed()

# Save the plot to the PDF file
dev.off()


```


```R
## Fig.3b

data <- read.table(file = "./dataset/LTRfamily_ERV1subfamily.txt",header=T,sep = '\t')

plot_frm <- data
plot_frm$Stage <- factor(plot_frm$Stage,levels = unique(plot_frm$Stage))
plot_frm$Subfamily <- factor(plot_frm$Subfamily, levels = c(
  "BHIKHARI_I", "BHIKHARI-5-I_DR", "BHIKHARI-3-I_DR", "ERV1-4_DR-I", "ERV1-3-I_DR", "ERV1_DR-I",
  "ERV1", "ERV1-N4-I_DR", "ERV1-N2-I_DR", "Gypsy","Pao","Copia", "Ngaro"

library(ggplot2)
library(tidyplot)
p <- ggplot(data = plot_frm, mapping = aes(x = Stage, y = Subfamily, size=TPM))+
  geom_point(shape=16,color='#f57c00')+
  scale_size_continuous(range = c(1,10))+
  facet_wrap(.~Type, nrow=2,scales = 'free_y')+
  theme_clean()
p

ggsave(filename = '../figure/bobule.pdf',plot = p,device = 'pdf',width = 5,height = 5.5,units = 'in')

## Fig. 3c 
heatmap_data <- read.table(file = "./dataset/bik_all_loci.log2TPM.txt",header=T, sep = '\t',row.names = 1)

blue_col <- colorRampPalette(c("#bdbdbd","#fafafa"))(5)
median_col <- colorRampPalette(c("#fafafa", "#fff3e0"))(60)
red_col <- colorRampPalette(c("#fff3e0","#e65100"))(35)
pal <- c(blue_col, red_col)

library(pheatmap)
p <- pheatmap(mat = heatmap_data, color = pal,scale='none',cluster_rows = T,cluster_cols = F, clustering_method = 'ward.D2', legend_breaks = seq(0,10,by=5), show_rownames = F,border_color = NA, treeheight_row = FALSE)
p

ggsave(filename = '../figure/heatmap.pdf',plot = p,device = 'pdf',width = 5,height = 3,units = 'in')
    
## Fig. 3d
library(data.table)
library(tidyverse)
input_tpm <- fread("./dataset/TE_trans.annotation.tpm.tsv",sep="\t",header=T,data.table = FALSE)
row.names(input_tpm) <- input_tpm$isoform
input_tpm_mat <- input_tpm[,-c(1:8)] %>% as.matrix()
input_tpm_mat <- input_tpm_mat[rowSums(input_tpm_mat)>0,]
library(pheatmap)
cluster_num <- 8
library(pheatmap)
p <- pheatmap(mat = input_tpm_mat, scale = "row", cutree_rows = cluster_num, show_rownames = F, cluster_cols = F, clustering_method = "ward.D",silent = T)
#------how to change cluster order
hclust.result<-p$tree_row
cluster = cutree(hclust.result,k=cluster_num)
cluster_name = paste("C",cluster,sep="")
names(cluster_name) = names(cluster)
hclust.result.order<-p$tree_row$order
cluster_name_ordered = cluster_name[hclust.result.order]
cluster_categroy = unique(cluster_name_ordered)
new_cluster_name = paste("cluster",1:length(cluster_categroy),sep="")
names(new_cluster_name) = cluster_categroy
anno_row = data.frame(cluster = cluster_name_ordered)
anno_row$cluster = new_cluster_name[as.character(anno_row$cluster)]
anno_row$cluster = factor(anno_row$cluster,
                          levels = paste("cluster",1:cluster_num,sep=""))

cluster_frm_tbl = data.frame(input_tpm[row.names(anno_row),],
                         cluster = anno_row[,"cluster"],check.names=F)
out_dir = "./dataset"
dir.create(out_dir,recursive = T)
write.table(cluster_frm_tbl,file.path(out_dir,"heatmap_tpm_clustered_by_isofrom.tsv"),sep="\t",quote=F,row.names=F)
row_gaps <- cumsum(table(cluster_frm_tbl$cluster))
input_tpm_ordered = input_tpm_mat[row.names(cluster_frm_tbl),]
library(RColorBrewer)
p = pheatmap(input_tpm_ordered, scale="row",cluster_row=T,cluster_col=F, clustering_method = "ward.D",
             show_rownames=F,annotation_row = anno_row, fontsize = 8, main="",cutree_rows = cluster_num,
             gaps_row = row_gaps,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

ggsave(filename = "heatmap_isoform_level.pdf",plot = p,device = "pdf",path = out_dir,width = 5,height = 5,units = "in")
    
## Fig. 3d cluster lines
    
library(reshape2)
cluster_frm <- reshape2::melt(data = cluster_frm_tbl[,c('isoform', 'cluster',colnames(input_tpm_ordered))],id.vars=c('isoform', 'cluster'),variable.name = "stage",value.name="TPM")
cluster_frm$stage <- factor(cluster_frm$stage, levels = colnames(input_tpm_mat))

cluster_median <- cluster_frm %>% group_by(cluster,stage) %>% summarise(median = median(log2(TPM+1)),sd = sd(log2(TPM+1)))
text_size=8
line_size=0.5/1.07
library(ggplot2)
l <- ggplot(data = cluster_median, mapping = aes(x = stage, y = median,group=cluster))+
#  geom_line(color="gray")+
  geom_line(color="red")+
#  geom_point(color="red",shape=16)+
  geom_ribbon(aes(ymin = median-sd, ymax = median+sd), alpha = 0.1)+
  facet_wrap(.~cluster,ncol=1)+
  theme_bw()
l
    ggsave(filename = "heatmap_isoform_cluster_line2.pdf",plot = l,device = "pdf",path = out_dir,width = 1.1,height = 7,units = "in")
    
## Fig. 3e cluster entropy
cluster_frm_subfamily <- table(cluster_frm_tbl$subfamily,cluster_frm_tbl$cluster)
cluster_frm_subfamily_num <- apply(cluster_frm_subfamily,1,function(x){
  length(x[which(x>0)])
})

cluster_frm_subfamily_cluster_tbl <- as.data.frame(cluster_frm_subfamily)
colnames(cluster_frm_subfamily_cluster_tbl) <- c('subfamily','cluster','subfam_isonum')

cluster_frm_subfamily_multi <- cluster_frm_subfamily[rowSums(cluster_frm_subfamily)>1,]
class(cluster_frm_subfamily_multi) <- "matrix"
cluster_frm_subfamily_entrop <- apply(cluster_frm_subfamily_multi,1,function(x){
  entropy::entropy(x/sum(x))})
cluster_frm_subfamily_multi_frm <- data.frame(subfamily = row.names(cluster_frm_subfamily_multi),
                                              entropy = cluster_frm_subfamily_entrop[row.names(cluster_frm_subfamily_multi)],
                                              cluster_frm_subfamily_multi)
write.table(cluster_frm_subfamily_multi_frm,file.path(out_dir,"subfamily_by_cluster_with_entropy.tsv"),sep="\t",quote=F,row.names=F)

#--visualization
plot_entropy <- cluster_frm_subfamily_entrop %>% as.data.frame()
colnames(plot_entropy) <- "entropy"

library(ggplot2)
library(tidyplot)
l <- ggplot(data = plot_entropy , mapping = aes(x = entropy))+
  geom_histogram(fill='#ab47bc',binwidth = 0.1,color="#000000")+
  theme_clean()+
  coord_flip()+
  scale_x_reverse()
l

ggsave(filename = "heatmap_isoform_cluster_entropy_bar.pdf", plot = l,device = "pdf",path = out_dir,width = 2,height = 2,units = "in")

cluster_frm_tbl2 <- left_join(x = cluster_frm_tbl,y = cluster_frm_subfamily_cluster_tbl, by=c('subfamily'='subfamily','cluster'='cluster'))
write.table(cluster_frm_tbl2,file.path(out_dir,"heatmap_tpm_clustered_by_isofrom_with_num.tsv"),sep="\t",quote=F,row.names=F)
```


```R
## Fig.5

## Fig.5a
library(ggplot2)
library(dplyr)

pdf("allstage_totalTPM_SR_LR.barplot.pdf",w=10,h=6)
data<-read.table("allstage_totalTPM_SR_LR.txt", header = TRUE, sep = "\t")
coeff <- 0.5
level_order <- factor(data$Stage, level = c("zygote","1-cell","2-cell","64-cell","128-cell","1k-cell","oblong","sphere","dome","shield", "75epiboly","1-4-somites","14-19-somites","20-25-somites","prim5","prim15","prim25","Long-pec","Protruding-mouth","day4","day5"))

ggplot(data, aes(x=level_order)) + 
  geom_bar( aes(y=SR / coeff), stat="identity", size=.1, fill="grey", color="grey", alpha=.9) +
  geom_line(aes(x=level_order,y=LR), size=1, color="#C00000",group = 1) + geom_point(aes(x=level_order,y=LR), color="#C00000",size=2)+
  #geom_line(aes(x=level_order, y=TPM / coeff), size=1, color="#B3697A",group = 1) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Long Read",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Short Read"))+ 
    theme_classic()+

  theme(
    axis.text.x = element_text(color = "black", size=16, angle=90), axis.text.y = element_text(color = "#D12E43", size=16), axis.text.y.right = element_text(color = "grey", size=16),axis.title.y = element_text(color = "#C00000", size=16),
    axis.title.y.right = element_text(color = "grey", size=16)
  )+

  ggtitle("Total TPM over early embryonic development")+theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5))

 # p+scale_x_discrete(limits=c("fertilized","cell1","cell64","cell1k","high","oblong","sphere","dome","30epiboly","50epiboly","shield"))
  
dev.off()





```
