library(tximportData)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(clusterProfiler)
library(DEGreport)
library(org.Mm.eg.db)
library(DOSE)
library(pathview)
library(tximport)
library(AnnotationDbi)
library(AnnotationHub)
library(ensembldb)
library(BiocFileCache)
library(vidger)

ah <- AnnotationHub()
yes
mus_ens <- query(ah, c("Mus musculus", "EnsDb"))
mus_ens
mus_ens <- mus_ens[["AH78811"]]
genes(mus_ens, return.type = "data.frame") %>% View()

# Create a transcript dataframe
txdb <- transcripts(mus_ens, return.type = "data.frame") %>%  dplyr::select(tx_id, gene_id)

txdbb <- txdbb[grep("ENST", txdbb$tx_id),]

# Create a gene-level dataframe
genedb <- genes(mus_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, symbol)

# Merge the two dataframes together
annotations <- inner_join(txdb, genedb)
annotations %>% View()

#give the location of samples

samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)

samples
dir <- system.file("extdata/salmon/salmon", package="tximportData")

files <- file.path(dir, samples$run, "quant.sf")
file.exists(files)
names(files) <- paste0("sample", 1:12)
files

# Run tximport
txii <- tximport(files, type="salmon", tx2gene=annotations[,c("tx_id", "gene_id")], countsFromAbundance="lengthScaledTPM", ignoreTxVersion = TRUE)

attributes(txii)

Treatment <- factor(c(rep("PBS", 3), rep("LPS_15min", 3), rep ("LPS_30min", 3), rep("LPS_4hrs", 3)))

Treatment <- factor(Treatment,levels = c("PBS", "LPS_15min", "LPS_30min", "LPS_4hrs"))

metaa <- data.frame(Treatment, row.names = colnames(txii$counts))


all(colnames(txii$counts) %in% rownames(metaa))
all(colnames(txii$counts) == rownames(metaa))

ddss <- DESeqDataSetFromTximport(txii, colData = metaa, design = ~Treatment)
ddss
View(counts(ddss))

ddss <- estimateSizeFactors(ddss)


ddss <- DESeq(ddss, test = "LRT", reduced = ~1)

normalized_countss <- counts(ddss, normalized=TRUE)

write.table(normalized_countss, file="data/normalized_countss.txt", sep="\t", quote=F, col.names=NA)

#-----PCA Plot

rowRanges(ddss)
nrow(ddss)

lambda <- 10^seq(from = -1, to = -2, length = 1000)

vsd <- vst(ddss, blind = FALSE)

head(assay(vsd), 3)

colData(vsd)

plotPCA(vsd, intgroup = "Treatment")

pcaData <- plotPCA(vsd, intgroup = "Treatment", returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x=PC1, y=PC2, color = Treatment)) + 
  geom_point(size = 6) + geom_point(alpha = 1) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_bw() + scale_colour_manual(values = c("#FF8C00", "#00ff00", "#FF00FF", "#CCCC00")) + theme(plot.title = element_text(size = 14), legend.title=element_text(size=18), legend.text=element_text(size=18), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) + ggtitle("")


#-----Scatterplot

vsScatterMatrix(data = ddss, d.factor = 'Treatment', type = 'deseq',  comp = NULL, title = TRUE, grid = TRUE, man.title = NULL)


#----Defining contrasts for each timepoint vs PBS----


contrast_4hrs <- c("Treatment", "LPS_4hrs", "PBS")
contrast_30min <- c("Treatment", "LPS_30min", "PBS")
contrast_15min <- c("Treatment", "LPS_15min", "PBS")


res_table4hrs <- results(ddss, contrast=contrast_4hrs, test = "Wald", alpha = 0.05, lfcThreshold = 0.58)

res_table30min <- results(ddss, contrast=contrast_30min, test = "Wald", alpha = 0.05, lfcThreshold = 0.58)

res_table15min <- results(ddss, contrast=contrast_15min, test = "Wald", alpha = 0.05, lfcThreshold = 0.58)


res_table4hrs <- as.data.frame(res_table4hrs)

res_table4hrs$significant <- ifelse(res_table4hrs$padj < .05 & abs(res_table4hrs$log2FoldChange) > 0.58, "Significant", NA)

ggplot(res_table4hrs, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=scales::squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="black", size=0.7, linetype="longdash") + labs(x="mean of normalized counts", y="log2 fold change") + scale_colour_manual(name="p.adj", values=("Significant"="red"), na.value="grey50", labels=c("Significant (4843)", "Not significant")) + theme_bw() + ggtitle("LPS 4hrs vs PBS - Cerebral Vessels") + theme(plot.title = element_text(size = 14), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_hline(yintercept=c(0.58, -0.58), linetype="dashed", color = "blue", size=0.5) + theme(legend.position="bottom")

res_table30min <- as.data.frame(res_table30min)

res_table30min$significant <- ifelse(res_table30min$padj < .05 & abs(res_table30min$log2FoldChange) > 0.58, "Significant", NA)

ggplot(res_table30min, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=scales::squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="black", size=0.7, linetype="longdash") + labs(x="mean of normalized counts", y="log2 fold change") + scale_colour_manual(name="p.adj", values=("Significant"="red"), na.value="grey50", labels=c("Significant (578)", "Not significant")) + theme_bw() + ggtitle("LPS 30min vs PBS - Cerebral Vessels") + theme(plot.title = element_text(size = 14), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_hline(yintercept=c(0.58, -0.58), linetype="dashed", color = "blue", size=0.5) + theme(legend.position="bottom")


res_table15min <- as.data.frame(res_table15min)

res_table15min$significant <- ifelse(res_table15min$padj < .05 & abs(res_table15min$log2FoldChange) > 0.58, "Significant", NA)

ggplot(res_table15min, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=scales::squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="black", size=0.7, linetype="longdash") + labs(x="mean of normalized counts", y="log2 fold change") + scale_colour_manual(name="p.adj", values=("Significant"="red"), na.value="grey50", labels=c("Significant (457)", "Not significant")) + theme_bw() + ggtitle("LPS 15min vs PBS - Cerebral Vessels") + theme(plot.title = element_text(size = 14), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_hline(yintercept=c(0.58, -0.58), linetype="dashed", color = "blue", size=0.5) + theme(legend.position="bottom")


padj.cutoff <- 0.05

res_table4hrs_tb <- res_table4hrs %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_table30min_tb <- res_table30min %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_table15min_tb <- res_table15min %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()


#----annotations

grch99annot <- annotations %>% 
  dplyr::select(gene_id, symbol) %>% 
  dplyr::distinct()

normalized_countss <- counts(ddss, normalized=T) %>%   data.frame() %>%
  rownames_to_column(var="gene") 


normalized_countss <- merge(normalized_countss, grch99annot, by.x="gene", by.y="gene_id")

normalized_countss <- normalized_countss %>%
  as_tibble()


sig4hrs <- res_table4hrs_tb %>%
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > 0.58)

sig4hrs_Txt <- merge(sig4hrs, grch99annot, by.x="gene", by.y="gene_id")

sig30min <- res_table30min_tb %>%
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > 0.58)

sig30min_Txt <- merge(res_table30min_tb, grch99annot, by.x="gene", by.y="gene_id")

sig15min <- res_table15min_tb %>%
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > 0.58)

sig15min_Txt <- merge(sig15min, grch99annot, by.x="gene", by.y="gene_id")


nrow(sig4hrs)
nrow(sig30min)
nrow(sig15min)

#For heatmaps, filtering genes by category

resp_to_lps_15min <- dplyr::filter(normalized_countss, symbol %in% c("Tnfaip3", "Nfkbia", "Tnf", "Il1b", "Ppbp", "Cxcl1", "Ptgs2", "Cxcl10", "Ccl12", "Gch1", "Cnr1", "Zfp36", "Jun", "Junb", "Cxcl2", "Jund"))


resp_to_lps_15min_Txt <- dplyr::filter(sig15min_Txt, symbol %in% c("Tnfaip3", "Nfkbia", "Tnf", "Il1b", "Ppbp", "Cxcl1", "Ptgs2", "Cxcl10", "Ccl12", "Gch1", "Cnr1", "Zfp36", "Jun", "Junb", "Cxcl2", "Jund"))


write.csv(resp_to_lps_15min_Txt, "data/resp_to_lps_15min_Txt.csv")


vessel15_vessel_morpho <- dplyr::filter(normalized_countss, symbol %in% c("Klf4", "Rgcc", "Adamts1", "Xdh", "Creb3l1", "Angpt4", "Sema4a", "Adgrb2", "Adgrb1", "Cxcl10", "Tgfb2", "Abcc8", "Klf2"))

vessel15_vessel_morpho_Txt <- dplyr::filter(sig15min_Txt, symbol %in% c("Klf4", "Rgcc", "Adamts1", "Xdh", "Creb3l1", "Angpt4", "Sema4a", "Adgrb2", "Adgrb1", "Cxcl10", "Tgfb2", "Abcc8", "Klf2"))


write.csv(vessel15_vessel_morpho_Txt, "data/vessel15_vessel_morpho_Txt.csv")


vessel30_vessel_dia <- dplyr::filter(normalized_countss, symbol %in% c("Alox12", "Tnf", "Ptgs2", "Manf", "P2ry2", "Dusp5", "Icam1", "Gch1", "Kcna5", "Per2", "Kdr", "Gpx1"))

vessel30_vessel_dia_Txt <- dplyr::filter(sig30min_Txt, symbol %in% c("Alox12", "Tnf", "Ptgs2", "Manf", "P2ry2", "Dusp5", "Icam1", "Gch1", "Kcna5", "Per2", "Kdr", "Gpx1"))


write.csv(vessel30_vessel_dia_Txt, "data/vessel30_vessel_dia_Txt.csv")

resp_to_lps_30min <- dplyr::filter(normalized_countss, symbol %in% c("Il10", "Tnfaip3", "Nfkbia", "Acod1", "Litaf", "Noct", "Tnf", "Rela", "Il6", "Il1rn", "Il1b", "Bcl10", "Cxcl1", "Ptgs2", "Trib1", "Nlrp3", "Cxcl10", "Ccl12", "Ccl2", "Sbno2", "Serpine1", "Arid5a", "Gch1", "Ripk2", "Zc3h12a", "Zfp36", "Cd14", "Cx3cr1", "Jun", "Junb", "Nod2", "Cebpb", "Cxcl2", "Jund"))


nfkb_sig_30min <- dplyr::filter(normalized_countss, symbol %in% c("Tnf", "Rela", "Il1b", "Tlr2", "Bcl10", "Nlrp3", "Amh", "Icam1", "Ripk2", "Prkd2", "Cx3cr1", "Nod2"))


tnf_prodn_30min <- dplyr::filter(normalized_countss, symbol %in% c("Ccl3", "Hspb1", "Il10", "Ccl4", "Tnfaip3", "Hspd1", "Tlr2", "Errfi1", "Ccl2", "Arid5a", "Akap12", "Thbs1", "Ripk2", "Zc3h12a", "Zfp36", "Trex1", "Cd14", "Cx3cr1", "Bcl3", "Nod2"))


cytokine_prodn_30min <- dplyr::filter(normalized_countss, symbol %in% c("Icosl", "Ccl3", "Hspb1", "Mcoln2", "Slamf6", "Il10", "Irf1", "Ccl4", "Rel", "Rgcc", "Tnf", "Rela", "Il6", "Hspd1", "Il1r1", "Il1rn", "Il1b", "Il1a", "Tlr2", "Bcl10", "Nr4a3", "Clec4e", "Ptgs2", "Nlrp3", "P2ry2", "Ccl2", "Tnfsf9", "Serpine1", "Arid5a", "Ddx60", "Egr1", "Akap12", "Thbs1", "Ripk2", "Prkd2", "Atf4", "Trim16", "C5ar1", "Cd14", "Bcl3", "Nod2", "Cebpg", "Cebpb", "Osm", "Gbp5"))


go_all_four_30min <- dplyr::filter(normalized_countss, symbol %in% c("Il10", "Tnfaip3", "Nfkbia", "Acod1", "Litaf", "Noct", "Tnf", "Rela", "Il6", "Il1rn", "Il1b", "Bcl10", "Cxcl1", "Ptgs2", "Trib1", "Nlrp3", "Cxcl10", "Ccl12", "Ccl2", "Sbno2", "Serpine1", "Arid5a", "Gch1", "Ripk2", "Zc3h12a", "Zfp36", "Cd14", "Cx3cr1", "Jun", "Junb", "Nod2", "Cebpb", "Cxcl2", "Jund", "Tlr2","Amh", "Icam1", "Prkd2", "Ccl3", "Hspb1", "Ccl4", "Hspd1", "Errfi1", "Akap12", "Thbs1", "Trex1", "Bcl3", "Icosl", "Mcoln2", "Slamf6", "Irf1", "Rel", "Rgcc", "Il1r1", "Il1a", "Nr4a3", "Clec4e", "P2ry2", "Tnfsf9", "Ddx60", "Egr1", "Atf4", "Trim16", "C5ar1", "Cebpg", "Osm", "Gbp5"))

go_all_four_30min_Txt <- dplyr::filter(sig30min_Txt, symbol %in% c("Il10", "Tnfaip3", "Nfkbia", "Acod1", "Litaf", "Noct", "Tnf", "Rela", "Il6", "Il1rn", "Il1b", "Bcl10", "Cxcl1", "Ptgs2", "Trib1", "Nlrp3", "Cxcl10", "Ccl12", "Ccl2", "Sbno2", "Serpine1", "Arid5a", "Gch1", "Ripk2", "Zc3h12a", "Zfp36", "Cd14", "Cx3cr1", "Jun", "Junb", "Nod2", "Cebpb", "Cxcl2", "Jund", "Tlr2","Amh", "Icam1", "Prkd2", "Ccl3", "Hspb1", "Ccl4", "Hspd1", "Errfi1", "Akap12", "Thbs1", "Trex1", "Bcl3", "Icosl", "Mcoln2", "Slamf6", "Irf1", "Rel", "Rgcc", "Il1r1", "Il1a", "Nr4a3", "Clec4e", "P2ry2", "Tnfsf9", "Ddx60", "Egr1", "Atf4", "Trim16", "C5ar1", "Cebpg", "Osm", "Gbp5"))


write.csv(go_all_four_30min_Txt, "data/go_all_four_30min_Txt.csv")


est_end_barrier_4hrs <- dplyr::filter(normalized_countss, symbol %in% c("Rock2", "Akap11", "Cldn1", "Tnf", "Il1b", "Abcb1b", "Add1", "Tnfrsf1a", "Msn", "Cdh5", "Myd88", "Icam1", "Ppp1r16b", "F11r", "Cldn5", "Sox18", "Rap2c", "Ezr", "Rap1b", "Rapgef2", "S1pr3", "Afdn", "Pde2a"))

est_end_barrier_4hrs_Txt <- dplyr::filter(sig4hrs_Txt, symbol %in% c("Rock2", "Akap11", "Cldn1", "Tnf", "Il1b", "Abcb1b", "Add1", "Tnfrsf1a", "Msn", "Cdh5", "Myd88", "Icam1", "Ppp1r16b", "F11r", "Cldn5", "Sox18", "Rap2c", "Ezr", "Rap1b", "Rapgef2", "S1pr3", "Afdn", "Pde2a"))

write.csv(est_end_barrier_4hrs_Txt, "data/est_end_barrier_4hrs_Txt.csv")


IEGs <- dplyr::filter(normalized_countss, symbol %in% c("Egr1", "Nr4a1", "Ier3", "Dusp1", "Egr3", "Fos", "Jun", "Fosb", "Junb", "Ier2", "Atf3", "Hspa1a")) 

IEGs_15min_Txt <- dplyr::filter(sig15min_Txt, symbol %in% c("Egr1", "Nr4a1", "Ier3", "Dusp1", "Egr3", "Fos", "Jun", "Fosb", "Junb", "Ier2", "Atf3", "Hspa1a"))

write.csv(IEGs_15min_Txt, "data/IEGs_15min_Txt.csv")

IEGs_30min_Txt <- dplyr::filter(sig30min_Txt, symbol %in% c("Egr1", "Nr4a1", "Ier3", "Dusp1", "Egr3", "Fos", "Jun", "Fosb", "Junb", "Ier2", "Atf3", "Hspa1a"))

write.csv(IEGs_30min_Txt, "data/IEGs_30min_Txt.csv")

IEGs_4hrs_Txt <- dplyr::filter(sig4hrs_Txt, symbol %in% c("Egr1", "Nr4a1", "Ier3", "Dusp1", "Egr3", "Fos", "Jun", "Fosb", "Junb", "Ier2", "Atf3", "Hspa1a"))

write.csv(IEGs_4hrs_Txt, "data/IEGs_4hrs_Txt.csv")


resp_to_lps_15min <- as.data.frame(resp_to_lps_15min)
row.names(resp_to_lps_15min) <- resp_to_lps_15min$symbol

go_all_four_30min <- as.data.frame(go_all_four_30min)
row.names(go_all_four_30min) <- go_all_four_30min$symbol

nrow(go_all_four_30min)

IEGs <- as.data.frame(IEGs)
row.names(IEGs) <- IEGs$symbol

ann_colors15 = list(
  Treatment = c(PBS = "#FF8C00", LPS_15min = "#00ff00"))

pheatmap(resp_to_lps_15min[2:7], 
         color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         annotation = (metaa),
         legend_breaks = c(-1.5,0,1.5),
         annotation_colors = ann_colors15,
         annotation_names_col = F,
         border_color = NA, 
         fontsize = 14, 
         scale = "row", 
         fontsize_row = 14, 
         height = 20)

annorow30 <- read.csv(file.path("annorow30min.csv"), header = TRUE, sep = ",", row.names = 1)

annorow30<-as.data.frame(annorow30)
annorow30$GO.0032496 <- factor(annorow30$GO.0032496, exclude = "")
annorow30$GO.0051092 <- factor(annorow30$GO.0051092, exclude = "")
annorow30$GO.0001819 <- factor(annorow30$GO.0001819, exclude = "")
annorow30$GO.0032640 <- factor(annorow30$GO.0032640, exclude = "")

ann_colors30 = list(
  Treatment = c(PBS = "#FF8C00", LPS_15min = "#00ff00", LPS_30min = "#FF00FF"), GO.0032496 = c("Response to lipopolysaccharide" = "#008080"), GO.0051092 = c("Positive regulation of NF-kappaB transcription factor activity" = "#F5DEB3"), GO.0001819 = c("Positive regulation of cytokine production" = "#A52A2A"), GO.0032640 = c("Tumor necrosis factor production" = "#FF7F84"))

pheatmap(go_all_four_30min[2:10], 
         color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         treeheight_row = 10,
         annotation = (metaa),
         annotation_names_col = F,
         legend_breaks = c(-2,0,2),
         annotation_row = annorow30,
         annotation_names_row = F,
         annotation_colors = ann_colors30,
         border_color = NA, 
         fontsize = 14, 
         scale = "row", 
         fontsize_row = 14, 
         height = 20,
         cellwidth = 50)



ann_colorsIEG = list(Treatment = c(PBS = "#FF8C00", LPS_15min = "#00ff00", LPS_30min = "#FF00FF", LPS_4hrs = "#CCCC00"))

pheatmap(IEGs[2:13], 
         color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         legend_breaks = c(-1.5,0,1.5),
         annotation = (metaa),
         annotation_names_col = F,
         annotation_colors = ann_colorsIEG,
         border_color = NA, 
         fontsize = 14, 
         scale = "row", 
         fontsize_row = 14, 
         height = 20)


vessel15_vessel_morpho <- as.data.frame(vessel15_vessel_morpho)
row.names(vessel15_vessel_morpho) <- vessel15_vessel_morpho$symbol

ann_colors15 = list(
  Treatment = c(PBS = "#FF8C00", LPS_15min = "#00ff00"))

pheatmap(vessel15_vessel_morpho[2:7], 
         color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         legend_breaks = c(-1,0,1),
         annotation = (metaa),
         annotation_names_col = F,
         annotation_colors = ann_colors15,
         border_color = NA, 
         fontsize = 14, 
         scale = "row", 
         fontsize_row = 14, 
         height = 20)


vessel30_vessel_dia <- as.data.frame(vessel30_vessel_dia)
row.names(vessel30_vessel_dia) <- vessel30_vessel_dia$symbol

ann_colors30 = list(
  Treatment = c(PBS = "#FF8C00", LPS_15min = "#00ff00", LPS_30min = "#FF00FF"))

pheatmap(vessel30_vessel_dia[2:10], 
         color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         legend_breaks = c(-1,0,1),
         annotation = (metaa),
         annotation_names_col = F,
         annotation_colors = ann_colors30,
         border_color = NA, 
         fontsize = 14, 
         scale = "row", 
         fontsize_row = 14, 
         height = 20)


est_end_barrier_4hrs <- as.data.frame(est_end_barrier_4hrs)
row.names(est_end_barrier_4hrs) <- est_end_barrier_4hrs$symbol

pheatmap(est_end_barrier_4hrs[2:13], 
         color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         legend_breaks = c(-1.5,0,1.5),
         annotation = (metaa),
         annotation_names_col = F,
         annotation_colors = ann_colorsIEG,
         border_color = NA, 
         fontsize = 14, 
         scale = "row", 
         fontsize_row = 14, 
         height = 20)


#---_GO enrichment using clusterprofiler

#pathway

res_15minensids <- inner_join(res_table15min_tb, annotations, by=c("gene"="gene_id"))    

res_30minensids <- inner_join(res_table30min_tb, annotations, by=c("gene"="gene_id"))

res_4hrsensids <- inner_join(res_table4hrs_tb, annotations, by=c("gene"="gene_id"))


## Create background dataset for hypergeometric testing using all genes tested for significance in the results 

all15mingenes <- as.character(res_15minensids$gene)
all30mingenes <- as.character(res_30minensids$gene)
all4hrsgenes <- as.character(res_4hrsensids$gene)

head(all15mingenes)

## Extract significant results

sig15 <- dplyr::filter(res_15minensids, padj < 0.05 & abs(log2FoldChange) > 0.58)
sig30 <- dplyr::filter(res_30minensids, padj < 0.05 & abs(log2FoldChange) > 0.58)
sig4 <- dplyr::filter(res_4hrsensids, padj < 0.05 & abs(log2FoldChange) > 0.58)


sig15_genes <- as.character(sig15$gene)
sig30_genes <- as.character(sig30$gene)
sig4hr_genes <- as.character(sig4$gene)


library(org.Mm.eg.db)



## Run GO enrichment analysis 
ego15 <- enrichGO(gene = sig15_genes, 
                  universe = all15mingenes,
                  keyType = "ENSEMBL",
                  OrgDb = org.Mm.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)

ego_top_15min <- ego15

ego_top_15min <- dplyr::filter(ego_top_15min, ID %in% c("GO:0002237", "GO:0032496", "GO:0050727", "GO:0002224", "GO:0019221", "GO:0031663", "GO:0032640", "GO:0002221", "GO:0007249", "GO:0051092"))



ego_top_15min

dotplot(ego_top_15min, showCategory=10)
emapplot(ego_top_15min)


ego30 <- enrichGO(gene = sig30_genes, 
                  universe = all30mingenes,
                  keyType = "ENSEMBL",
                  OrgDb = org.Mm.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)

ego_top_30min <- ego30

ego_top_30min  

ego_top_30min <- dplyr::filter(ego_top_30min, ID %in% c("GO:0002237", "GO:0032496", "GO:0050727", "GO:0002224", "GO:0019221", "GO:0031663", "GO:0032640", "GO:0002221", "GO:0007249", "GO:0051092"))

dotplot(ego_top_30min, showCategory = 10)           

emapplot(ego_top_30min)
emapplot(ego30)


ego4hr <- enrichGO(gene = sig4hr_genes, 
                   universe = all4hrsgenes,
                   keyType = "ENSEMBL",
                   OrgDb = org.Mm.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)


ego_top_4hr <- ego4hr

ego_top_4hr <- dplyr::filter(ego_top_4hr, ID %in% c("GO:0002237", "GO:0032496", "GO:0050727", "GO:0002224", "GO:0019221", "GO:0031663", "GO:0032640", "GO:0002221", "GO:0007249", "GO:0051092"))

dotplot(ego_top_4hr, showCategory=10)
emapplot(ego_top_4hr)



ego_vessel_related_4hr <- dplyr::filter(ego4hr, ID %in% c("GO:2001233", "GO:0034330", "GO:0061028", "GO:0007044", "GO:1901888", "GO:0034332", "GO:0007045", "GO:0090109", "GO:1901889", "GO:2000351", "GO:0043297", "GO:0120192", "GO:1903392", "GO:0016264", "GO:2000353"))

dotplot(ego_vessel_related_4hr, showCategory=10)

emapplot(ego_vessel_related_4hr)

## Output results from GO analysis to a table
cluster_summary15 <- data.frame(ego15)

write.csv(cluster_summary15, "data/clusterProfiler_15.csv")

cluster_summary30 <- data.frame(ego30)

write.csv(cluster_summary30, "data/clusterProfiler_30.csv")

cluster_summary4hr <- data.frame(ego4hr)

write.csv(cluster_summary4hr, "data/clusterProfiler_4hr.csv")

