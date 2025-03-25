#Create dotplot for the Fst results 
#Y-axis will be the -log p value where we plot the largest p-value 
#X-axis will be the virus names 
#Size of the dot will correspond to how many population pairs the virus was associated with 
library(stringr)
library(Biostrings)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(xlsx)
neg_controls <- c('GRANTHAM', 'DOOLITTLE', 'SURFACE')
output_dir <- '~/Figure_4/Data/Fst_analysis/Results' #directory where all Fst comparison results are stored 
plot_dir <- '~/Figure_4/Plots/Fst_analysis'
dir.create(plot_dir, recursive = T)
pop_pairs <- str_remove(list.dirs(output_dir, recursive = F), paste0(output_dir,'/'))
cut_off <- 75
analysis_df <- data.frame('Virus' = numeric(), 'Organism' = character(), 'Significance' = character(),
                            'log_p_val' = numeric(), "Pop_pair" = character())
for (p in 1:length(pop_pairs)) {
    pop_df <- fread(paste0(output_dir, '/', pop_pairs[p], '/', cut_off, '/', 'INTERFACE_p_values.txt'))
    pop_df$Pop_pair <- pop_pairs[p]
    analysis_df <- rbind(analysis_df, pop_df[, c('Virus', 'organism', "Significance", "log_p_val", "Pop_pair")])
}
result_df <- analysis_df %>%
  filter(Significance == '***' & log_p_val >= 1.3) %>%  # Filter only rows where Significance is '***'
  group_by(Virus, organism) %>%
    ungroup()
#get the viruses that were significant in at least one population pair 
sig_viruses <- unique(result_df$Virus)

plot_df <- analysis_df[analysis_df$Virus %in% sig_viruses, ]
# Create the dot plot
 # Create a new column for colors based on significance and population pair
 plot_df$log_p_val <- ifelse( plot_df$log_p_val >= 5, 5,  plot_df$log_p_val)
 
 plot_df <- plot_df %>%
   mutate(color = ifelse(log_p_val >= 1.3, as.character(Pop_pair), "grey"))
 
 #write.xlsx(unique(plot_df[, c('Virus', 'organism')]), file = '/ix/djishnu/Priyamvada/viral-human-ppi/viral-human-ppi-R-scripts/publication/Fst_final/plot_df.xlsx', row.names = F)
 labels <- read.xlsx('/ix/djishnu/Priyamvada/viral-human-ppi/viral-human-ppi-R-scripts/publication/Fst_final/plot_df_labels_final.xlsx', sheetIndex = 1)
 labels <- labels[labels$Keep == 'Yes',]
 label_vector <- setNames(str_replace_all(labels$Label, fixed(" "), ""), labels$Virus)
 order_of_labels <-  read.xlsx('/ix/djishnu/Priyamvada/viral-human-ppi/viral-human-ppi-R-scripts/publication/Fst_final/plot_df_labels_final.xlsx', sheetIndex = 3)
 order_of_labels$Labels <- str_replace_all(order_of_labels$Label, fixed(" "), "")
 
 plot_df <- plot_df[plot_df$Virus %in% labels$Virus,]
 # Define a color palette for your population pairs
 color_palette <- c("AFRvEAS" = '#ff595e', 
                    "AFRvSAS" = '#ff924c', 
                    "EASvSAS" = '#ffca3a', 
                    "EURvAFR" = '#8ac926', 
                    "EURvEAS" = '#1982c4', 
                    "EURvSAS" = '#6a4c93', 
                    "grey" = "grey")
int_plot <- ggplot(
  plot_df, 
  aes(x = as.character(Virus), y = log_p_val, shape = Pop_pair)) +  
  geom_point(aes(color = color), size = 4) +
  guides(color = 'none') +
  scale_color_manual(values = color_palette) + # Ensure correct color mapping
  scale_shape_manual(name = "Population Pair",
                     labels = c("African-East Asian", "African-South Asian", "East Asian-South Asian",
                                'European-African', 'European-East Asian', 'European-South Asian'),
                     values = c(16, 17, 15, 3, 7, 8)) +
  labs(
    x = "Virus",
    y = "Negative Log of p-value"
  ) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size = 0.75), 
    plot.title = element_text(face = 'bold'), 
    strip.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(face = 'bold', size = 14, color = 'black'), 
    strip.background = element_blank(),
    strip.text.x = element_text(size = 14, color = 'black'),
    axis.text.x = element_text(size = 12, color = 'black', angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 12, color = 'black')
  ) + 
  scale_x_discrete(labels = label_vector, limits = factor(order_of_labels$Virus)) +
  geom_hline(yintercept=1.3, linetype="dashed", 
             color = "red", size=0.5)
#dev.off()
ggsave(filename = 'INTERFACE_plots.pdf', plot = int_plot, device = 'pdf', path = plot_dir, width = 10, height = 5, units = 'in')

#generate plots for negative controls 
neg_analysis_list <- list()
for (a in 1:length(neg_controls)) {
  neg_analysis_df <- data.frame('Virus' = numeric(), 'Organism' = character(), 'Significance' = character(),
                                'log_p_val' = numeric(), "Pop_pair" = character())
  for (p in 1:length(pop_pairs)) {
    pop_df <- fread(paste0(output_dir, '/', pop_pairs[p], '/', cut_off, '/', neg_controls[a], '_p_values.txt'))
    pop_df$Pop_pair <- pop_pairs[p]
    neg_analysis_df <- rbind(neg_analysis_df, pop_df[, c('Virus', 'organism', "log_p_val", "Pop_pair")])
  }
  neg_analysis_df <- neg_analysis_df[neg_analysis_df$Virus %in% labels$Virus,]
  neg_analysis_df <- neg_analysis_df %>%
                      mutate(color = ifelse(log_p_val >= 1.3, as.character(Pop_pair), "grey"))
  neg_plot <- ggplot(
    neg_analysis_df, 
    aes(x = as.character(Virus), y = log_p_val, shape = Pop_pair)) +  
    geom_point(aes(color = color), size = 4) +
    guides(color = 'none') +
    scale_color_manual(values = color_palette) + # Ensure correct color mapping
    scale_shape_manual(name = "Population Pair",
                       labels = c("African-East Asian", "African-South Asian", "East Asian-South Asian",
                                  'European-African', 'European-East Asian', 'European-South Asian'),
                       values = c(16, 17, 15, 3, 7, 8)) +
    labs(
      x = "Virus",
      y = "Negative Log of p-value"
    ) +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black", size = 0.75), 
      plot.title = element_text(face = 'bold'), 
      strip.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(face = 'bold', size = 12, color = 'black'), 
      strip.background = element_blank(),
      strip.text.x = element_text(size = 14, color = 'black'),
      axis.text.x = element_text(size = 12, color = 'black', angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 12, color = 'black')
    ) + 
    scale_x_discrete(labels = label_vector, limits = factor(order_of_labels$Virus)) +
    ylim(0, 5) +
    geom_hline(yintercept=1.3, linetype="dashed", 
               color = "red", size=0.5)
  ggsave(filename = paste0(neg_controls[a], '_plots.pdf'), plot = neg_plot, device = 'pdf', path = plot_dir, width = 10, height = 5, units = 'in')
} 

