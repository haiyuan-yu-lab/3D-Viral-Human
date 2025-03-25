library(data.table)
library(stringr)
library(ggpubr)
library(rstatix)
library(dplyr)
plot_df <- read.table('//ix/djishnu/Priyamvada/viral-human-ppi/Github/3D-Viral-Human/Figure_5/Data/Vir_Hu_Hu/Mutation_rate_comparison.txt', header = T)
output_dir <- '/ix/djishnu/Priyamvada/viral-human-ppi/Github/3D-Viral-Human/Figure_5/Plots/Vir_Hu_Hu'
stat_test <- plot_df %>%
  wilcox_test(Dice_index ~ Virus, p.adjust.method = 'BH') %>%
  add_significance() %>%
  add_xy_position(x = "Virus" , dodge = 0.8 , scales = "free")
median_jaccard_df <- plot_df %>%
  group_by(Virus) %>%
  summarize(median_jaccard = median(Dice_index, na.rm = TRUE))
# Calculate the IQR, lower, and upper bounds
Q1 <- quantile(plot_df$Dice_index, 0.25)
Q3 <- quantile(plot_df$Dice_index, 0.75)
IQR_value <- IQR(plot_df$Dice_index)

# Define lower and upper bounds
lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value

# Filter the dataframe to remove outliers
plot_df_filtered <- plot_df %>%
  filter(Dice_index >= lower_bound & Dice_index <= upper_bound)
colors <- c('#00429d','#93003a')
names(colors) <- c('>=1e-04', '<1e-04')
bxp <- ggplot(
  plot_df_filtered, 
  aes(x = Virus, y = Dice_index)) +
  geom_violin(trim=T, aes(fill = Virus)) +
  geom_boxplot(width=0.075, notch = T, outliers = F) +
  scale_fill_manual(values = colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.75), 
        plot.title = element_text(face = 'bold'), axis.text = element_text(size = 14, color = 'black'), 
        axis.title = element_text(face = 'bold', size = 16, color = 'black'), legend.position = 'none', strip.background = element_blank(),
        strip.text.x = element_text(size = 14, color = 'black'), ) +
  labs(x = "", y = "EXTENT OF SHAREDNESS") 
pdf(paste0(output_dir, '/', 'Mutation_comparison_publication.pdf'), height = 5, width = 6)
print(bxp + 
        stat_pvalue_manual(
          stat_test, y.position = 0.55, tip.length = 0, size = 6,
          label = "{p.signif}"
        ))
dev.off()
