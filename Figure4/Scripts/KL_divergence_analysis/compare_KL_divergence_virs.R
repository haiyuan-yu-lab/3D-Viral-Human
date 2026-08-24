#For each virus, compare KL divergence values for interface vs non-interface residues.
library(tidyverse)
library(ggridges)
library(rstatix)
library(data.table)
library(stringr)

score_cutoff <- 0.3
PPV_cutoff   <- 0.6

parent_dir <- "~/3D-Viral-Human/Figure4/"
data_dir   <- file.path(parent_dir, "Data/KL_divergence_analysis")
source(file.path(data_dir, "Scripts/KL_divergence_analysis/functions_KL_divergence.R"))

pred_df <- fread(file.path('/path/to/3Dinteractome.csv'))  # Replace with the actual path to interactome data 
pred_df <- pred_df[pred_df$score >= score_cutoff, ]
pred_df <- pred_df[pred_df$PPV   >= PPV_cutoff,   ]

taxa     <- c(385599, 11676, 211044, 2697049)
names(taxa) <- c("H3N2", "HIV", "H1N1", "COVID")
virnames <- c("H3N2", "HIV", "H1N1", "COVID")

kl_annot_all <- data.frame("Vir_Prot" = character(), "Gene" = character(),
                           "POS" = integer(), "KL" = numeric(),
                           "Interface" = character(), "Vir" = character())

for (t in seq_along(taxa)) {
  tax_id    <- taxa[t]
  vir_prots <- pred_df[pred_df$Viral_taxonomy_id == tax_id, ]
  
  kl_dir   <- file.path(data_dir, names(taxa)[t], "KL_vals/")
  kl_files <- list.files(kl_dir)
  splits   <- str_split(kl_files, "_")
  protein_uids <- vapply(splits, `[`, character(1), 2)
  protein_uids <- str_remove(protein_uids, "\\.txt$")
  gene_names   <- vapply(splits, `[`, character(1), 1)
  
  gene_map <- data.frame("Viral_protein_uid" = protein_uids,
                         "Viral_gene_name"   = gene_names,
                         stringsAsFactors = FALSE)
  
  kl_annot_df <- data.frame("Vir_Prot" = character(), "Gene" = character(),
                            "POS" = integer(), "KL" = numeric(),
                            "Interface" = character())
  
  for (p in seq_along(kl_files)) {
    prot <- str_trim(gene_map[p, 1])
    gene <- gene_map[p, 2]
    print(paste0("Processing ", gene))
    
    prot_interface <- vir_prots[vir_prots$Viral_protein_uid == prot, ]
    keep <- !is.na(prot_interface$Viral_ires_plddt50) &
      prot_interface$Viral_ires_plddt50 != "-"
    prot_interface <- prot_interface[keep, , drop = FALSE]
    
    if (nrow(prot_interface) == 0) {
      message("  [WARN] no interface annotations for protein: ", prot)
      next
    }
    
    pldtt_ires <- extract_positions(prot_interface$Viral_ires_plddt50)
    
    kl_df <- read.table(paste0(kl_dir, "/", gene, "_", prot, ".txt"), header = TRUE)
    kl_df <- kl_df %>%
      mutate(Interface = case_when(POS %in% pldtt_ires ~ "I", TRUE ~ "NI")) %>%
      filter(Interface %in% c("I", "NI"))
    
    if (!(all(kl_df$Interface == "I") || all(kl_df$Interface == "NI"))) {
      kl_annot_df <- rbind(kl_annot_df,
                           data.frame("Vir_Prot" = rep(prot, nrow(kl_df)),
                                      "Prot"      = rep(prot, nrow(kl_df)),
                                      "Gene"      = rep(gene, nrow(kl_df)),
                                      "POS"       = kl_df$POS,
                                      "KL"        = kl_df$KL,
                                      "Interface" = kl_df$Interface))
    } else {
      print(paste0("Not sufficient observations for not interface in ", gene))
    }
  }
  # ── combine dataframes ───────────────────────────────────────────────────────────────────
  comp_levels <- c("I", "NI")
  kl_annot_df <- kl_annot_df %>%
    mutate(Interface = factor(Interface, levels = comp_levels))
  kl_annot_df$Vir <- names(taxa)[t]
  kl_annot_all    <- rbind(kl_annot_all, kl_annot_df)
}

# ── all-virus stats ────────────────────────────────────────────────────────────
comp_levels  <- c("I", "NI")
kl_annot_all <- kl_annot_all %>%
  mutate(Interface = factor(Interface, levels = comp_levels))

stat_test_all <- kl_annot_all %>%
  group_by(Vir) %>%
  wilcox_test(KL ~ Interface, p.adjust.method = "BH") %>%
  add_significance() %>%
  add_xy_position(x = "Interface", dodge = 0.8, scales = "free")
print(paste0("P value for wilcox test is ", stat_test_all$p))

whisker_df_all <- kl_annot_all %>%
  group_by(Vir, Interface) %>%
  summarise(Q1 = quantile(KL, 0.25, na.rm = TRUE),
            Q3 = quantile(KL, 0.75, na.rm = TRUE),
            IQR = Q3 - Q1,
            upper_whisker = max(KL[KL <= (Q3 + 1.5 * IQR)], na.rm = TRUE),
            .groups = "drop") %>%
  group_by(Vir) %>%
  summarise(y_position = max(upper_whisker, na.rm = TRUE) * 0.9)

stat_test_all <- stat_test_all %>% left_join(whisker_df_all, by = "Vir")

# Cliff's delta + median difference per virus
cliff_df_all <- kl_annot_all %>%
  group_by(Vir) %>%
  cliff_delta(KL ~ Interface) %>%
  ungroup() %>%
  mutate(delta = round(delta, 2),
         cliff_label = paste0("\u03B4=", delta)) %>%
  left_join(whisker_df_all, by = "Vir") 

print(cliff_df_all %>% select(Vir, cliff_label))

# Merge into stat_test for bracket annotation
stat_test_all <- stat_test_all %>%
  left_join(cliff_df_all %>% select(Vir, cliff_label),
            by = "Vir") %>%
  mutate(final_label = paste0(p.signif, " (", cliff_label, ")"))

# ── all-virus boxplot ──────────────────────────────────────────────────────────
bxp_all <- ggplot(kl_annot_all, aes(x = Interface, y = KL)) +
  geom_boxplot(aes(fill = Interface)) +
  scale_fill_manual(values = c("#8856a7", "#9ebcda")) +
  labs(x = "RESIDUE TYPE", y = "KL DIVERGENCE VALUES") +
  facet_wrap(~Vir, scales = "free_x") +
  theme(panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        panel.background  = element_blank(),
        axis.line         = element_line(colour = "black", size = 0.75),
        plot.title        = element_text(face = "bold"),
        axis.text         = element_text(size = 10, color = "black"),
        axis.title        = element_text(face = "bold", size = 12, color = "black"),
        legend.position   = "none",
        strip.background  = element_blank(),
        strip.text.x      = element_text(size = 10, color = "black"))

out_dir <- file.path(data_dir, "Plots")
dir.create(out_dir, recursive = TRUE)

# ── ridge plots ────────────────────────────────────────────────────────────────
theme_set(theme_void())
theme_flip <- theme(
  axis.text.x        = element_text(color = "black", face = "plain", size = 10,
                                     margin = margin(t = 4)),
  axis.text.y        = element_text(color = "black", face = "plain", size = 10,
                                     hjust = 1, margin = margin(b = 75)),
  axis.line.x        = element_line(color = "black", size = 0.5),
  panel.grid.major.x = element_line(color = "grey90", size = .1),
  panel.grid.major.y = element_blank(),
  legend.position    = "top",
  legend.text        = element_text(size = 16),
  plot.title         = element_text(face = "bold", size = 14, margin = margin(t = 5, l = 60, b = 5)),
  axis.title.y       = element_text(margin = margin(l = 20, t = 10, b = 5), size = 10),
  axis.title.x       = element_text(margin = margin(l = 20, t = 10, b = 5), size = 10),
  plot.background    = element_rect(fill = "white", color = "white"),
  plot.margin        = margin(t = 20, r = 20, b = 20, l = 20)
)


my_pal <- c('I' = '#4b1b5a', 'NI' = '#FCF6F5FF')
label_vector <- c("NI" = "Not Interface Residues", "I" = "Interface Residues")
vir_labels   <- c("COVID" = "SARS-CoV-2", "H1N1" = "H1N1", "H3N2" = "H3N2", "HIV" = "HIV")

# Build per-virus facet label: "SARS-CoV-2\nδ=X  ∆Mdn=Y"
ridge_annot <- cliff_df_all %>%
  select(Vir, cliff_label) %>%
  deframe()

vir_labels_with_effect <- setNames(
  paste0(vir_labels, "\n", ridge_annot[names(vir_labels)]),
  names(vir_labels)
)

kl_annot_all_plot <- kl_annot_all[kl_annot_all$KL <= 2, ]

y <- ggplot(kl_annot_all_plot,
            aes(KL, fct_rev(Interface), color = Interface, fill = Interface)) +
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = c(.07, .07)) +
  scale_color_manual(values = my_pal, guide = "none") +
  scale_fill_manual(values = my_pal, guide = "none") +
  theme_flip

ridge_plot <- y +
  ggridges::stat_density_ridges(
    quantile_lines = TRUE, quantiles = 2,
    color = "black", alpha = .8, size = 1.0
  ) +
  labs(y = NULL, x = "KL Divergence Values",
       title = "Human-Viral Protein Interactions") +
  scale_y_discrete(expand = c(.1, .1), labels = label_vector) +
  facet_wrap(~Vir, scales = "free_x",
             labeller = as_labeller(vir_labels_with_effect))

print(ridge_plot)
ggsave(filename = "All_vir_ridgeplot.pdf", plot = ridge_plot,
       device = cairo_pdf, path = out_dir, height = 7.5, width = 10, units = "in")