library(ggplot2)
library(ggsignif)
library(patchwork)
library(scales)

##################################

# pLDDT distribution plots
plddts <- read.csv('C:/Users/crtuser/Documents/CrAss PHD/Crassvirales_ICTV_proposal_2024/Structural_Taxonomy/data/pLDDTs_proteins.csv')

plddts$protein <- gsub(".pdb$", "", plddts$protein)  # Remove '.pdb'

df <- plddts

# Histogram of pLDDT scores
hist_plot <- ggplot(df, aes(x = pLDDT)) +
  geom_histogram(bins = 10, fill = "steelblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "pLDDT Scores Distribution", x = "pLDDT Score", y = "Frequency")
hist_plot

# Scatter plot of pLDDT vs protein length
protein_len_plddt_plot <- ggplot(df, aes(x = protein_length, y = pLDDT)) +
  geom_point(color = "steelblue", size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "pLDDT Score vs Protein Length", x = "Protein Length", y = "pLDDT Score")
protein_len_plddt_plot


# Genome length plots
genome_lengths <- read.csv('C:/Users/crtuser/Documents/CrAss PHD/Crassvirales_ICTV_proposal_2024/Structural_Taxonomy/data/4083_total_folded_per_genome.csv')
genome_lengths <- genome_lengths[order(genome_lengths$proteins_folded), ]
genome_lengths$genome_length <- as.numeric(genome_lengths$genome_length)


genome_lengths[!complete.cases(genome_lengths[, c("genome_length", "proteins_folded")]), ]


genomes_prot_folded_plot <- ggplot(genome_lengths, aes(x = genome_length, y = proteins_folded)) +
  geom_point(color = "steelblue", size = 1.5, alpha = 0.7) +
  theme_minimal() +
  scale_x_continuous(labels = scientific) +
  labs(
    title = "Total Proteins Folded by Genome Length",
    x = "Genome Length",
    y = "Proteins Folded"
  )

genomes_prot_folded_plot


combined_plots <- (
  protein_len_plddt_plot | hist_plot | genomes_prot_folded_plot
) 

library(cowplot)
final_plot <- ggdraw(combined_plots) +
  draw_plot_label(
    label = c("A", "B", "C"),
    x = c(0.015, 0.35, 0.67),  # x positions for each subplot
    y = c(0.99, 0.99, 0.99),   # y positions (1 = top)
    size = 10
  )

final_plot

ggsave("C:/Users/crtuser/Documents/CrAss PHD/Crassvirales_ICTV_proposal_2024/Structural_Taxonomy/plots/pdbs_macro_metadata.png", final_plot, width = 10, height = 5, dpi = 300)










