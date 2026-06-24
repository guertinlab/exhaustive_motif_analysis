bar=read.table('peak.with.motif.sum.txt', sep = "", header=F)
colnames(bar)=c("enriched_motifs", "peak_numbers")
bar$dum.x="peak number"
bar$enriched_motifs<- factor(bar$enriched_motif, levels = c("other", "motif_11_161bp", "motif_10_161bp",  "motif_9_161bp", "motif_8_161bp", "motif_7_161bp",  "motif_6_161bp", "motif_5_161bp", "motif_2_161bp", "motif_1_161bp"))

#coherence check
sum(bar$peak_numbers) 

library("lattice") 
library(RColorBrewer)

pdf('peak_number_with_or_without_ZNF143_motif_variant.pdf', width=5, height=8)


par.settings = list(superpose.line = list(col=c("red", "black"), lwd=2), strip.background=list(col="grey85"))

my.settings <- list(
  #superpose.polygon=list(col=c(colorRampPalette(c("red","pink"))(6),colorRampPalette(c("blue","light blue"))(3), "light grey"), border="transparent"),
  #superpose.polygon=list(col=c("red","orange","yellow","lightgreen","skyblue","blue","#FF00FF","brown","#CC6600", "#006666", "lightgrey"), border="transparent"),
  superpose.polygon=list(col=c("light grey", brewer.pal(11, "Spectral"), "red"),  border="transparent"),
  strip.background=list(col="grey80", cex = 0.6),
  strip.border=list(col="black")
)
print(barchart(peak_numbers ~ dum.x,         
               data = bar,
               groups = enriched_motifs,
               stack = TRUE,
               auto.key = list(space = "right", lines=F, points=F, rectangles = TRUE, cex = 1, reverse.rows = TRUE),
               #auto.key=list(space="right"),
               #scales = list(x = list(rot = 45)),
               ylab = "number of ChIP-seq peaks",
               xlab = "",
               par.settings = my.settings,
               panel=function(x, y, ...) {
                 panel.barchart(x, y, ...) 
                 # Add text labels
                 # panel.text(x, cumsum(y)-y, labels = y, pos = 3, cex = 0.6)
               }
)
)
dev.off()