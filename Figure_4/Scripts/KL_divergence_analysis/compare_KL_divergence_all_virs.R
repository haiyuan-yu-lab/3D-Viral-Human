#for each protein compare the interface and not interface positions 
library(tidyverse)    
library(data.table)   
library(stringr)      
library(rstatix)      
library(ggridges)     

# Load viral-human interaction data and filter by confidence score >= 0.10
all_pred_files <- fread('~/Figure_4/Data/KL_divergence_analysis/viral-human_predictions.txt')
all_pred_files <- all_pred_files[all_pred_files$score >= 0.10, ]


# Define taxonomy IDs for viruses of interest
taxa <- c(385599, 11676, 211044, 2697049)
names(taxa) <- c('H3N2', 'HIV', 'H1N1', 'COVID')

# Initialize empty dataframes to store KL annotation and p-values
KL_annot_all <- data.frame('Vir_Prot' = character(), 'Gene' = character(), 'POS' = integer(), 'KL' = numeric(), 'Interface' = character(), 'Vir' = character())
KL_p_vals_df <- data.frame('Vir' = character(), 'P_value' = numeric())

# Loop over each virus
for (t in 1:length(taxa)) {
  print(names(taxa)[t])
  tax_id <- taxa[t]
  
  # Filter interactions for the current virus
  H1N1_prots <- all_pred_files[all_pred_files$Viral_taxonomy_id == tax_id, ]
  
  # Define directory where KL divergence values are stored
  KL_file_dir <- paste0('~/Figure_4/Data/KL_divergence_analysis/', names(taxa)[t], '/KL_vals/')
  KL_files <- list.files(KL_file_dir)
  
  # Extract gene and protein names from filenames
  H1N1_gene_map <- data.frame(
    'Viral_protein_uid' = str_remove(sapply(str_split(KL_files, '_'), function(x) x[2]), '.txt'), 
    'Viral_gene_name' = sapply(str_split(KL_files, '_'), function(x) x[1])
  )
  
  # Initialize empty dataframe to store KL values for current virus
  KL_annnot_df <- data.frame('Vir_Prot' = character(), 'Gene' = character(), 'POS' = integer(), 'KL' = numeric(), 'Interface' = character())
  
  # Loop over each protein for the virus
  for (p in 1:length(KL_files)) {
    prot <- str_trim(H1N1_gene_map[p, 1])
    gene <- H1N1_gene_map[p, 2]
    
    # Extract interface annotations for the protein
    prot_interface <- H1N1_prots[H1N1_prots$Viral_protein_uid == prot, ]
    prot_interface <- prot_interface[!Viral_ires_plddt50 == "-", ]
    prot_interface <- prot_interface[!Viral_ires == "-", ]
    
    # Extract high-confidence interface and all interface residues
    pldtt_ires <- unique(as.numeric(unlist(strsplit(prot_interface$Viral_ires_plddt50, ','))))
    all_ires <- unique(as.numeric(unlist(strsplit(prot_interface$Viral_ires, ','))))
    low_conf_ires <- all_ires[!all_ires %in% pldtt_ires]
    
    # Read KL divergence values
    KL_df <- read.table(paste0(KL_file_dir, KL_files[which(KL_files == paste0(gene, '_', prot, '.txt'))]), header = TRUE)
    
    # Annotate each position as Interface (I) or Non-Interface (NI)
    KL_df <- KL_df %>% mutate('Interface' = case_when(POS %in% pldtt_ires ~ 'I', TRUE ~ 'NI'))
    KL_df <- KL_df[KL_df$Interface == 'I' | KL_df$Interface == 'NI', ]
    
    # Include only proteins that have both I and NI positions
    if (!(all(KL_df$Interface == 'I') | all(KL_df$Interface == 'NI'))) {
      KL_annnot_df <- rbind(KL_annnot_df, data.frame('Vir_Prot' = rep(prot, nrow(KL_df)), 'Prot' = rep(prot, nrow(KL_df)), 'Gene' = rep(gene, nrow(KL_df)),
                                                     'POS' = KL_df$POS, 'KL' = KL_df$KL, 'Interface' = KL_df$Interface))
    } else {
      print(paste0('Not sufficient observations for not interface in ', gene))
      KL_annnot_df <- rbind(KL_annnot_df, data.frame('Vir_Prot' = rep(prot, nrow(KL_df)), 'Prot' = rep(prot, nrow(KL_df)), 'Gene' = rep(gene, nrow(KL_df)),
                                                     'POS' = KL_df$POS, 'KL' = KL_df$KL, 'Interface' = KL_df$Interface))
    }
  }
  
  # Add virus label and merge with global dataframe
  KL_annnot_df$Vir <- names(taxa)[t]
  KL_annot_all <- rbind(KL_annot_all, KL_annnot_df)
}

# Perform Wilcoxon test comparing KL divergence between Interface vs Non-Interface for each virus
stat.test <- KL_annot_all %>%
  group_by(Vir) %>%
  wilcox_test(KL ~ Interface, p.adjust.method = "BH") %>%
  add_significance() %>%
  add_xy_position(x = "Vir", dodge = 0.8, scales = "free")

print(paste0('P value for wilcox test is ', stat.test$p))

# Define custom themes for plotting
theme_set(theme_void())
theme_update(
  axis.text.x = element_text(color = "black", size = 10, margin = margin(t = 4)),
  axis.text.y = element_text(color = "black", size = 10, hjust = 1, margin = margin(r = 6)),
  axis.line.x = element_line(color = "black", size = 0.5),
  plot.title = element_text(face = 'bold', size = 10, margin = margin(l = 30)),
  panel.grid.major.y = element_line(color = "grey90", size = .2),
  plot.background = element_rect(fill = "white", color = "white"),
  axis.title.y = element_text(margin = margin(t = 10), size = 10)
)

theme_flip <- theme(
  axis.text.x = element_text(face = "plain", size = 10),
  axis.text.y = element_text(face = "plain", size = 10, margin = margin(b = 75)),
  panel.grid.major.x = element_line(color = "grey90", size = .1),
  panel.grid.major.y = element_blank(),
  legend.position = "top", 
  legend.text = element_text(size = 16),
  plot.title = element_text(face = 'bold', size = 14, margin = margin(t = 5, l = 60, b = 5)),
  axis.title.y = element_text(margin = margin(l = 20, t = 10, b = 5), size = 10),
  axis.title.x = element_text(margin = margin(l = 20, t = 10, b = 5), size = 10),
  plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
)

# Set color palette and labels
my_pal <- c('I' = '#4b1b5a', 'NI' = '#FCF6F5FF')
label_vector <- c('NI' = 'Not Interface Residues', 'I' = 'Interface Residues')

# Filter extreme KL values
KL_annot_all_plot <- KL_annot_all[KL_annot_all$KL <= 2, ]
# Base ggplot
y <- ggplot(KL_annot_all_plot, aes(KL, fct_rev(Interface), color = Interface, fill = Interface)) + 
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = c(.07, .07)) +
  scale_color_manual(values = my_pal, guide = "none") +
  scale_fill_manual(values = my_pal, guide = "none") +
  theme_flip

# Set virus labels for plot facets
vir_labels <- c('COVID' = 'SARS-CoV-2', 'H1N1' = 'H1N1', 'H3N2' = 'H3N2', 'HIV' = 'HIV')

# Generate ridge plot
ridge_plot <- y + 
  ggridges::stat_density_ridges(
    quantile_lines = TRUE, quantiles = 2, 
    color = "black", alpha = .8, size = 1.0
  ) + 
  labs(y = NULL, x = "KL Divergence Values", title = 'Human-Viral Protein Interactions') +
  scale_y_discrete(expand = c(.1, .1), labels = label_vector) + 
  facet_wrap(~factor(Vir, levels = c('COVID','H1N1','HIV', 'H3N2')), scales = 'free_x', labeller = as_labeller(vir_labels))

# Display plot
print(ridge_plot)

# Save plot as PDF
ggsave(filename = 'All_viruses_compare_rearranged.pdf', plot = ridge_plot, 
       device = 'pdf', path = '~/Figure_4/Plots/KL_divergene_analysis', height = 7.5, width = 10, units = 'in')
