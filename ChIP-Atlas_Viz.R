library(ggplot2)
library(cowplot)
library(reshape2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

#Run the script on terminal with bed file as 1st arg and peak file as 2nd arg
args <- commandArgs(trailingOnly=TRUE)

peakCount <- nrow(read.table(args[2], header = FALSE, sep = "\t"))
prefix <- gsub(".bed", "", args[1])
overlaps <- read.table(args[1], sep = "\t", header = FALSE)
countsTable <- as.data.frame(xtabs(~V7, overlaps))
colnames(countsTable)[1] <- "ID"
countsTable$Percent <- round((countsTable$Freq/as.numeric(peakCount))*100,0)
filtered_ids <- countsTable %>% filter(countsTable$Percent >= 20) %>% droplevels()
id_factor <- data.frame( overlaps$V7, overlaps$V8, overlaps$V9)
colnames(id_factor) <- c("ID", "Factor", "Cell.Group")
merged_table <- merge(filtered_ids, id_factor, by = "ID", all.x = TRUE)
merged_table <- merged_table[!duplicated(merged_table),]
rownames(merged_table) <- NULL
merged_table$Cell.Type <- gsub(".*\\(@\\s*(.*?)\\s*\\).*", "\\1", merged_table$Factor)
merged_table$Factor <- gsub(" \\(.*", "", gsub("Name=", "",merged_table$Factor))
merged_table$Cell.Group <- gsub("Cell group=", "", merged_table$Cell.Group)
merged_table <- merged_table %>% filter(!grepl("Epitope tags", Factor))
merged_table <- merged_table %>% filter(!grepl("Biotin", Factor))
write.table(merged_table, file = paste0(prefix, "_table.tsv"), sep = "\t", row.names = FALSE)
num_colors <- length(unique(merged_table$Cell.Group))
color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(num_colors)
ggplot(data = merged_table, mapping = aes(x = Cell.Type, y = Factor, color = Cell.Group, size = Percent)) + geom_jitter(position = position_jitter(width = 0, height = 0)) + theme_bw() +theme(legend.text = element_text(size=7, color="black"), axis.title = element_text(size=10,face="bold"), legend.title = element_text(size = 7, color = "black"), axis.text = element_text(color = "black",face="bold", size = 10), axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+labs(title="Gain in B-ALL",y = "Factors", x = "Cell Type", color = "Cell Group") + scale_size_continuous(breaks = c(10,40,70,100), limits = c(0,100)) + scale_color_manual(values = c("Blood" = "red", "Embryo" = "cyan3", "Adipocyte" = "violet", "Breast" = "pink", "Others" = "green4")) + guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(paste0(prefix, "_plot.pdf"), height = 15, width = 20, units = "cm")

