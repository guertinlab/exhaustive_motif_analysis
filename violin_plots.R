library(ggplot2)


motif1 <- read.table("ZNF143_meme_classic_PSWM_in_peaks_1.bed")
df1 <- data.frame(motif="motif_1", intensity = motif1$V5)

motif2 <- read.table("ZNF143_meme_classic_PSWM_in_peaks_2.bed")
df2 <- data.frame(motif="motif_2", intensity = motif2$V5)

motif3 <- read.table("ZNF143_meme_de_PSWM_in_peaks_5.bed")
df3 <- data.frame(motif="motif_3", intensity = motif3$V5)

motif4 <- read.table("ZNF143_meme_de_PSWM_in_peaks_6.bed")
df4 <- data.frame(motif="motif_4", intensity = motif4$V5)

motif5 <- read.table("ZNF143_meme_de_PSWM_in_peaks_7.bed")
df5 <- data.frame(motif="motif_5", intensity = motif5$V5)

motif6 <- read.table("ZNF143_meme_de_PSWM_in_peaks_8.bed")
df6 <- data.frame(motif="motif_6", intensity = motif6$V5)

motif7 <- read.table("ZNF143_meme_de_PSWM_in_peaks_9.bed")
df7 <- data.frame(motif="motif_7", intensity = motif7$V5)

motif8 <- read.table("ZNF143_meme_de_PSWM_in_peaks_10.bed")
df8 <- data.frame(motif="motif_8", intensity = motif8$V5)

motif9 <- read.table("ZNF143_meme_de_PSWM_in_peaks_11.bed")
df9 <- data.frame(motif="motif_9", intensity = motif9$V5)

df <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)

ggplot(df, aes (x=motif, y=intensity)) + geom_violin() + labs(title="Peak Subset Strength", x = "ZNF143 Motifs", y="ChIP Peak Intensity") + geom_boxplot()
