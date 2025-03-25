library(tidyverse)
library(stringr)
library(ggpubr)
library(rstatix)

jaccard_df <- read.table('/ix/djishnu/Priyamvada/viral-human-ppi/Github/3D-Viral-Human/Figure_5/Data/Vir_Hu_Vir/HIV_HPV_comparison.txt',
                     header = T)
output_dir <- '/ix/djishnu/Priyamvada/viral-human-ppi/Github/3D-Viral-Human/Figure_5/Plots/Vir_Hu_Vir/'
colors <- c('#00429d','#ffffe0', '#93003a')
names(colors) <- c('HIV:HIV', 'HIV:High risk', 'High risk:High risk')

stat.test <- jaccard_df %>%
  wilcox_test(Dice_index ~ Comparison_type,  p.adjust.method = "BH") %>%
  add_significance() %>%
  add_xy_position(x = "Comparison_type" , dodge = 0.8 , scales = "free")
# Create a box plot
bxp <- ggplot(
  jaccard_df, 
  aes(x = Comparison_type, y = Dice_index)) +
  geom_violin(trim=T, aes(fill = Comparison_type, color = 'black')) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) + 
  scale_x_discrete(breaks=c('HIV:HIV', 'HIV:High risk', 'High risk:High risk'),
                   labels = c('HIV:HIV', 'HIV:High risk HPV', 'High risk HPV:High risk HPV')) +
  geom_boxplot(width=0.075, notch = T, outliers = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.title = element_text(size = 14, face = 'bold'), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11)) +
  labs(x = "", y = "EXTENT OF SHAREDNESS") 
stat.test$y.position <- c(0.95, 1.01, 0.95)
pdf(paste0(output_dir, '/', 'HIV_HPV_compare.pdf'), height = 3.9, width = 6)
print(bxp +  
        stat_pvalue_manual(
          stat.test, y.position = stat.test$y.position,
          label = "{p.adj.signif}", tip.length = 0, bracket.shorten = 0.10, size = 4, 
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) 
)
dev.off()