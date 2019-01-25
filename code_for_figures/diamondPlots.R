

#############  DIAMOND PLOTS
# In these plots, each gene is represented by a stack of diamonds corresponding to that gene’s associated accessible chromatin regions. The bottom-most peak in this stack corresponds to the log2 fold change in expression of the gene. The diamonds are colored according to the accessibility change of their ATAC-seq peak with blue indicating closing and red indicating opening. The color scale was based on the rank-order of the peak accessibility changes. In Extended Data Fig. 6d, the color scale ranges from a log2 fold change of −3.92 to 4.96 (L14/L7).


library(openxlsx)
library(data.table)
library(magrittr)
library(dplyr)

#all.atac.res.df <- read.table("all.atac.res.df.txt" , sep = "\t", header=T, fill=T, stringsAsFactors=F)
all.atac.res.df <- read.xlsx("all.atac.res.df.xlsx")
all.rna.res.df <- read.xlsx("/scratchLocal/paz2005/andrew/sept11/rna/all.rna.res.df.xlsx")
row.names(all.rna.res.df) <- all.rna.res.df[,1]


genes_of_interest <- c("Fas", "Il2", "Tnfrsf25", "Tnf", "Nr4a1")
WT_vs_KO_de_res <-read.table("../rna/deseq2_WT_vs_KO_shrink.txt", header=T, sep="\t")
WT_vs_KO_de_res <- data.frame(SYMBOL=row.names(WT_vs_KO_de_res), log2FoldChange=WT_vs_KO_de_res$log2FoldChange)
WT_vs_KO_de_res <- WT_vs_KO_de_res[WT_vs_KO_de_res$SYMBOL %in% genes_of_interest,] 
WT_vs_KO_de_with_atac <- merge(all.atac.res.df, WT_vs_KO_de_res, by="SYMBOL")
up_and_down_with_atac_sorted <- WT_vs_KO_de_with_atac[order(WT_vs_KO_de_with_atac$log2FoldChange), ]
up_and_down_with_atac_sorted$SYMBOL = factor(up_and_down_with_atac_sorted$SYMBOL, levels=c(unique( up_and_down_with_atac_sorted$SYMBOL)))
WT_vs_KO_de_top25_with_atac_df <- as.data.frame(up_and_down_with_atac_sorted %>% group_by(SYMBOL) %>% dplyr::mutate(id = seq(0,0.3, length.out=n() ) ))

pdf("A_diamondPlot_WT_vs_KO_genes_of_interest.pdf")
WT_vs_KO_de_top25_with_atac_df <- as.data.frame(up_and_down_with_atac_sorted %>% group_by(SYMBOL) %>% dplyr::mutate(id = seq(0,0.3, length.out=n() ) ))
ggplot(WT_vs_KO_de_top25_with_atac_df, aes(SYMBOL,as.numeric( log2FoldChange+id), color=WT_vs_KO.Fold, label=SYMBOL)) +
  geom_point(shape=18, size=2.5, alpha = .9) + theme_classic() +
  theme_classic(base_size=10) + xlab("") +  scale_colour_gradient2("ATAC log2FC", low="blue", mid="white", high="red") + ylab("expression log2 fold change (RNA-seq)")  + theme(axis.text.x = element_text(angle = 90, hjust = 1))  + ggtitle("WT vs. KO") + geom_hline(yintercept=0, linetype="dotted")
dev.off()


## WT_vs_KO
top_n=25
WT_vs_KO_de_res <-read.table("../rna/deseq2_WT_vs_KO_shrink.txt", header=T, sep="\t")
#WT_vs_KO_de_res <- subset(all.rna.res.df, WT_vs_KO.padj < 0.05)
WT_vs_KO_de_res <- data.frame(SYMBOL=row.names(WT_vs_KO_de_res), log2FoldChange=WT_vs_KO_de_res$log2FoldChange)
WT_vs_KO_de_res <- WT_vs_KO_de_res[which(!is.na(WT_vs_KO_de_res$log2FoldChange)),]
WT_vs_KO_de_res <- WT_vs_KO_de_res[order(WT_vs_KO_de_res$log2FoldChange), ]
WT_vs_KO_de_pos_logfc <- WT_vs_KO_de_res[WT_vs_KO_de_res$log2FoldChange > 0,]
WT_vs_KO_de_neg_logfc <- WT_vs_KO_de_res[WT_vs_KO_de_res$log2FoldChange < 0,]
WT_vs_KO_de_pos_logfc <- WT_vs_KO_de_pos_logfc[WT_vs_KO_de_pos_logfc$SYMBOL %in% all.atac.res.df$SYMBOL,]
WT_vs_KO_de_neg_logfc <- WT_vs_KO_de_neg_logfc[WT_vs_KO_de_neg_logfc$SYMBOL %in% all.atac.res.df$SYMBOL,]
WT_vs_KO_de_top25_up <- tail(WT_vs_KO_de_pos_logfc, top_n)
WT_vs_KO_de_top25_down <- head(WT_vs_KO_de_neg_logfc, top_n)
WT_vs_KO_de_top25_up_with_atac <- merge(all.atac.res.df, WT_vs_KO_de_top25_up, by="SYMBOL")
WT_vs_KO_de_top25_down_with_atac <- merge(all.atac.res.df, WT_vs_KO_de_top25_down, by="SYMBOL")

up_and_down_with_atac <- rbind(WT_vs_KO_de_top25_down_with_atac, WT_vs_KO_de_top25_up_with_atac)
up_and_down_with_atac <- up_and_down_with_atac[,c("SYMBOL", "log2FoldChange", "WT_vs_KO.Fold")]
up_and_down_with_atac_sorted <- up_and_down_with_atac[order(up_and_down_with_atac$log2FoldChange), ]
up_and_down_with_atac_sorted$SYMBOL = factor(up_and_down_with_atac_sorted$SYMBOL, levels=c(unique( up_and_down_with_atac_sorted$SYMBOL)))

pdf("A_diamondPlot_WT_vs_KO.pdf", width=10)
WT_vs_KO_de_top25_with_atac_df <- as.data.frame(up_and_down_with_atac_sorted %>% group_by(SYMBOL) %>% dplyr::mutate(id = seq(0,0.3, length.out=n() ) ))
ggplot(WT_vs_KO_de_top25_with_atac_df, aes(SYMBOL,as.numeric( log2FoldChange+id), color=WT_vs_KO.Fold, label=SYMBOL)) +
  geom_point(shape=18, size=2.5, alpha = .9) + theme_classic() +
  theme_classic(base_size=10) + xlab("") +  scale_colour_gradient2("ATAC log2FC", low="blue", mid="white", high="red") + ylab("expression log2 fold change (RNA-seq)")  + theme(axis.text.x = element_text(angle = 90, hjust = 1))  + ggtitle("WT vs. KO") + geom_hline(yintercept=0, linetype="dotted")
dev.off()




## SV40_vs_OT1
top_n=25
SV40_vs_OT1_de_res <-read.table("../rna/deseq2_SV40_vs_OT1_shrink.txt", header=T, sep="\t")
SV40_vs_OT1_de_res <- data.frame(SYMBOL=row.names(SV40_vs_OT1_de_res), log2FoldChange=SV40_vs_OT1_de_res$log2FoldChange)
SV40_vs_OT1_de_res <- SV40_vs_OT1_de_res[which(!is.na(SV40_vs_OT1_de_res$log2FoldChange)),]
SV40_vs_OT1_de_res <- SV40_vs_OT1_de_res[order(SV40_vs_OT1_de_res$log2FoldChange), ]
SV40_vs_OT1_de_pos_logfc <- SV40_vs_OT1_de_res[SV40_vs_OT1_de_res$log2FoldChange > 0,]
SV40_vs_OT1_de_neg_logfc <- SV40_vs_OT1_de_res[SV40_vs_OT1_de_res$log2FoldChange < 0,]
SV40_vs_OT1_de_pos_logfc <- SV40_vs_OT1_de_pos_logfc[SV40_vs_OT1_de_pos_logfc$SYMBOL %in% all.atac.res.df$SYMBOL,]
SV40_vs_OT1_de_neg_logfc <- SV40_vs_OT1_de_neg_logfc[SV40_vs_OT1_de_neg_logfc$SYMBOL %in% all.atac.res.df$SYMBOL,]
SV40_vs_OT1_de_top25_up <- tail(SV40_vs_OT1_de_pos_logfc, top_n)
SV40_vs_OT1_de_top25_down <- head(SV40_vs_OT1_de_neg_logfc, top_n)
SV40_vs_OT1_de_top25_up_with_atac <- merge(all.atac.res.df, SV40_vs_OT1_de_top25_up, by="SYMBOL")
SV40_vs_OT1_de_top25_down_with_atac <- merge(all.atac.res.df, SV40_vs_OT1_de_top25_down, by="SYMBOL")

up_and_down_with_atac <- rbind(SV40_vs_OT1_de_top25_down_with_atac, SV40_vs_OT1_de_top25_up_with_atac)
up_and_down_with_atac <- up_and_down_with_atac[,c("SYMBOL", "log2FoldChange", "SV40_vs_OT1.Fold")]
up_and_down_with_atac_sorted <- up_and_down_with_atac[order(up_and_down_with_atac$log2FoldChange), ]
up_and_down_with_atac_sorted$SYMBOL = factor(up_and_down_with_atac_sorted$SYMBOL, levels=c(unique( up_and_down_with_atac_sorted$SYMBOL)))

pdf("B_diamondPlot_SV40_vs_OT1.pdf", width=10)
SV40_vs_OT1_de_top25_with_atac_df <- as.data.frame(up_and_down_with_atac_sorted %>% group_by(SYMBOL) %>% dplyr::mutate(id = seq(0,0.3, length.out=n() ) ))
ggplot(SV40_vs_OT1_de_top25_with_atac_df, aes(SYMBOL,as.numeric( log2FoldChange+id), color=SV40_vs_OT1.Fold, label=SYMBOL)) +
  geom_point( shape=18, size=2.5, alpha = .9) + theme_classic() +
  theme_classic(base_size=12) + xlab("") +  scale_colour_gradient2("ATAC log2FC", low="blue", mid="white", high="red") + ylab("expression log2 fold change (RNA-seq)")  + theme(axis.text.x = element_text(angle = 90, hjust = 1))  + ggtitle("SV40 vs. OT1") + geom_hline(yintercept=0, linetype="dotted")
dev.off()


## KO_vs_OT1
top_n=25
KO_vs_OT1_de_res <-read.table("../rna/deseq2_KO_vs_OT1_shrink.txt", header=T, sep="\t")
KO_vs_OT1_de_res <- data.frame(SYMBOL=row.names(KO_vs_OT1_de_res), log2FoldChange=KO_vs_OT1_de_res$log2FoldChange)
KO_vs_OT1_de_res <- KO_vs_OT1_de_res[which(!is.na(KO_vs_OT1_de_res$log2FoldChange)),]
KO_vs_OT1_de_res <- KO_vs_OT1_de_res[order(KO_vs_OT1_de_res$log2FoldChange), ]
KO_vs_OT1_de_pos_logfc <- KO_vs_OT1_de_res[KO_vs_OT1_de_res$log2FoldChange > 0,]
KO_vs_OT1_de_neg_logfc <- KO_vs_OT1_de_res[KO_vs_OT1_de_res$log2FoldChange < 0,]
KO_vs_OT1_de_pos_logfc <- KO_vs_OT1_de_pos_logfc[KO_vs_OT1_de_pos_logfc$SYMBOL %in% all.atac.res.df$SYMBOL,]
KO_vs_OT1_de_neg_logfc <- KO_vs_OT1_de_neg_logfc[KO_vs_OT1_de_neg_logfc$SYMBOL %in% all.atac.res.df$SYMBOL,]
KO_vs_OT1_de_top25_up <- tail(KO_vs_OT1_de_pos_logfc, top_n)
KO_vs_OT1_de_top25_down <- head(KO_vs_OT1_de_neg_logfc, top_n)
KO_vs_OT1_de_top25_up_with_atac <- merge(all.atac.res.df, KO_vs_OT1_de_top25_up, by="SYMBOL")
KO_vs_OT1_de_top25_down_with_atac <- merge(all.atac.res.df, KO_vs_OT1_de_top25_down, by="SYMBOL")

up_and_down_with_atac <- rbind(KO_vs_OT1_de_top25_down_with_atac, KO_vs_OT1_de_top25_up_with_atac)
up_and_down_with_atac <- up_and_down_with_atac[,c("SYMBOL", "log2FoldChange", "KO_vs_OT1.Fold")]
up_and_down_with_atac_sorted <- up_and_down_with_atac[order(up_and_down_with_atac$log2FoldChange), ]
up_and_down_with_atac_sorted$SYMBOL = factor(up_and_down_with_atac_sorted$SYMBOL, levels=c(unique( up_and_down_with_atac_sorted$SYMBOL)))



pdf("C_diamondPlot_KO_vs_OT1.pdf", width=10)
KO_vs_OT1_de_top25_with_atac_df <- as.data.frame(up_and_down_with_atac_sorted %>% group_by(SYMBOL) %>% dplyr::mutate(id = seq(0,0.3, length.out=n() ) ))
ggplot(KO_vs_OT1_de_top25_with_atac_df, aes(SYMBOL,as.numeric( log2FoldChange+id), color=KO_vs_OT1.Fold, label=SYMBOL)) +
  geom_point( shape=18, size=2.5, alpha = .9) + theme_classic() +
  theme_classic(base_size=12) + xlab("") +  scale_colour_gradient2("ATAC log2FC", low="blue", mid="white", high="red") + ylab("expression log2 fold change (RNA-seq)")  + theme(axis.text.x = element_text(angle = 90, hjust = 1))  + ggtitle("KO vs. OT1") + geom_hline(yintercept=0, linetype="dotted")
dev.off()

