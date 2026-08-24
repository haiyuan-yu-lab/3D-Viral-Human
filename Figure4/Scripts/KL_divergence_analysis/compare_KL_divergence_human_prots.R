#plotting for paper using code from https://z3tt.github.io/beyond-bar-and-box-plots/#box-plot-x-violin-plot-x-beeswarm-plot
library(tidyverse)
library(colorspace)
library(ggridges)
library(ggpubr)
library(rstatix)
library(data.table)
library(stringr)
KL_all_df <- fread('/ix/djishnu/Priyamvada/viral-human-ppi/Github/3D-Viral-Human/Figure_4/Data/KL_divergence_analysis/Human_prots_KL_vals.txt')

theme_set(theme_void())

theme_update(
  axis.text.x = element_text(color = "black", size = 12, 
                             margin = margin(t = 4)),
  axis.text.y = element_text(color = "black", size = 12, hjust = 1, 
                             margin = margin(r = 6)),
  axis.line.x = element_line(color = "black", size = 0.5),
  plot.title = element_text(face = 'bold', size = 14, margin = margin(l = 30)),
  panel.grid.major.y = element_line(color = "grey90", size = .2),
  plot.background = element_rect(fill = "white", color = "white"),
  axis.title.y = element_text(margin = margin(t = 10),
                              size = 16)
)
theme_flip <-
  theme(
    axis.text.x = element_text(face = "plain", size = 12),
    axis.text.y = element_text(face = "plain", size = 12, margin = margin(b = 75)),
    panel.grid.major.x = element_line(color = "grey90", size = .1),
    panel.grid.major.y = element_blank(),
    legend.position = "top", 
    legend.text = element_text(size = 16),
    plot.title = element_text(face = 'bold', size = 14, margin = margin(t = 5, l = 60)),
    axis.title.y = element_text(margin = margin(l = 20, t = 10, b = 5),
                                size = 14),
    axis.title.x = element_text(margin = margin(l = 20, t = 10, b = 5),
                                size = 14), 
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20) 
  )

## custom colors
my_pal <- c('I' = '#8856a7', 'NI' = '#9ebcda')
label_vector <- c('NI' = 'Not Interface Residues', 'I' = 'Interface Residues')
stat.test <- KL_all_df %>%
  wilcox_test(KL ~ Interface) %>%
  add_significance() %>%
  add_xy_position(x = "Interface" , dodge = 0.8 , scales = "free")

  
 #ridge plot 
 y <- ggplot(KL_all_df, aes(KL, fct_rev(Interface), color = Interface, fill = Interface)) + 
   coord_cartesian(clip = "off") +
   scale_y_discrete(expand = c(.07, .07)) +
   scale_color_manual(values = my_pal, guide = "none") +
   scale_fill_manual(values = my_pal, guide = "none") +
   theme_flip 
   
 
ridge_plot <- y + 
   ggridges::stat_density_ridges(
     quantile_lines = TRUE, quantiles = 2, 
     color = "black", alpha = .8, size = 1.0
   ) + labs(y = NULL, x = "KL Divergence Values", title = 'Human Protein Interactions') +
  scale_y_discrete(expand = c(.1, .1), labels = label_vector) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) + 
  xlim(0, 2)


ggsave(filename = 'Human_I_v_NI.pdf', plot = ridge_plot, 
       device = 'pdf', path = '/ix/djishnu/Priyamvada/viral-human-ppi/Github/3D-Viral-Human/Figure_4/Plots/KL_divergene_analysis', height = 5.5, width = 7, units = 'in')
