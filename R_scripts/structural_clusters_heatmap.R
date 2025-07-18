library(readr)  
library(tidyr)  
library(ggplot2)
library(dplyr)
library(RColorBrewer)

setwd('C:/Users/crtuser/Documents/Structural_Taxonomy/data/')

heatmap_data <- read.csv('../most_conserved_clusters_by_order_percentage.csv')

heatmap_long <- pivot_longer(
  heatmap_data,
  cols = -cluster_label,
  names_to = "Order",
  values_to = "value"
)

desired_order <- c(
 'Phage major capsid protein E | YP_010082820.1',
 "Major virion structural protein Mu-1/Mu-1C (M2) | YP_008240652.1",
 'Terminase RNaseH-like domain | YP_009620569.1',
 'Terminase RNaseH-like domain | YP_009783023.1',
 
 'Phage tail tube protein | YP_008058468.1',
 'Phage tail tube protein | YP_009302436.1',
 'Tail tubular protein | YP_009098374.1',
 
 'Tail terminator (DUF3168) | YP_010063459.1', 
 'Tail terminator (DUF3168) | YP_003347506.1',
 'T4-like virus tail tube protein gp19 | YP_009147612.1',
 
 'Phage portal protein, SPP1 Gp6-like: | YP_009802963.1',
 'Phage Connector (GP10) | YP_009908791.1',
 'Phage portal protein | YP_009807042.1'
)

heatmap_long <- heatmap_long %>%
  filter(
    cluster_label %in% desired_order
  )

setdiff(unique(heatmap_long$cluster_label), desired_order)

setdiff(desired_order,unique(heatmap_long$cluster_label))

desired_order_prefixed <- paste0(seq_along(desired_order), ". ", desired_order)

desired_order <- rev(desired_order)
desired_order_prefixed <- rev(desired_order_prefixed)
name_map <- setNames(desired_order_prefixed, desired_order)

# Set factor with proper labels and ordering
heatmap_long$cluster_label <- factor(
  heatmap_long$cluster_label,
  levels = desired_order,
  labels = desired_order_prefixed  # This handles the label renaming
)


specific_label <- "Unclassified_Order"

# Reorder the 'Family_member' column so 'specific_label' is the last
heatmap_long$Order <- factor(heatmap_long$Order, 
                                     levels = c(setdiff(unique(heatmap_long$Order), specific_label), specific_label))


colors <- brewer.pal(n = 9, name = "Blues")
last_8_colors <- colors[2:9]
last_8_colors[1] <- "white"
color_palette <- last_8_colors


# Define breaks and labels for heatmap bins
value_breaks <- c(-Inf, 0, 0.01, 0.05, 0.15, 0.50, 0.85, 0.99, Inf)
labels <- c("0", "<1%", "1–5%", "5–15%", "15–50%", "50–85%", "85–99%", "Core (99%+)")

heatmap_long$value_category <- cut(heatmap_long$value, 
                                   breaks = value_breaks, 
                                   include.lowest = TRUE, 
                                   labels = labels)


y_positions <- c(11.5,9.5,6.5,3.5)  # Example positions

heatmap_long_clean <- na.omit(heatmap_long)

# Create the heatmap with defined breaks and categories
p <- ggplot(heatmap_long, aes(x = Order, y = cluster_label)) +
  geom_tile(aes(fill = value_category), color = "grey") +
  scale_y_discrete(labels = name_map) +  # Use name_map to add prefixes to y-axis labels
  scale_fill_manual(values = color_palette, na.value = "white") +  # Map categories to colors
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.4),  # Adjust vjust as needed
        axis.ticks.x = element_line()) +
  labs(title = "Most Conserved Structural Cluster Presence by Order",
       x = "Order",
       y = "Structural Cluster",
       fill = "% Genomes with\n Protein in Cluster") +
  # Add horizontal lines only between specified y-axis categories
  geom_hline(yintercept = y_positions, color = "black", linewidth = 0.5)+
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 18, hjust = 1),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))

p

ggsave("../plots/heatmap_most_conserved_structures.tiff", plot = p, width = 8, height = 8, dpi = 600)


############ Create same plot for PFAM abundance #########

library(tidyverse)
library(RColorBrewer)
library(tibble)

# === Load data ===
df <- read.csv('../data/embeddings_PFAMS_easy_cluster.csv')  # Should contain: Pfam1, Function, DE1, genome_accn, Order

# === Step 1: Create Pfam × Order matrix of genome counts ===
pfam_counts <- df %>%
  group_by(Order, Pfam1) %>%
  summarise(count = n_distinct(genome_accn), .groups = 'drop')

# === Step 2: Spread into wide format ===
pfam_wide <- pfam_counts %>%
  pivot_wider(names_from = Order, values_from = count, values_fill = 0)

# === Step 3: Calculate median presence per Pfam across all Orders ===
pfam_medians <- pfam_wide %>%
  rowwise() %>%
  mutate(median_presence = median(c_across(-Pfam1))) %>%
  ungroup()

# === Step 4: Add Function for grouping ===
pfam_medians <- pfam_medians %>%
  left_join(df %>% select(Pfam1, Function) %>% distinct(), by = "Pfam1")

# === Step 5: Select top 1 most shared Pfam per Function ===
top_pfams <- pfam_medians %>%
  group_by(Function) %>%
  slice_max(order_by = median_presence, n = 1, with_ties = FALSE) %>%
  ungroup()

selected_pfams <- top_pfams$Pfam1

# === Step 6 (Corrected): Normalize counts to % of genomes per order ===
order_totals <- df %>%
  group_by(Order) %>%
  summarise(total_genomes = n_distinct(genome_accn), .groups = "drop")

pfam_norm <- df %>%
  filter(Pfam1 %in% selected_pfams) %>%
  group_by(Order, Pfam1) %>%
  summarise(count = n_distinct(genome_accn), .groups = 'drop') %>%
  left_join(order_totals, by = "Order") %>%
  mutate(Percent = count / total_genomes) %>%
  select(Order, Pfam1, Percent) %>%
  pivot_wider(names_from = Pfam1, values_from = Percent, values_fill = 0)



# === Step 7: Create labels: Function | Pfam | DE1 ===
label_df <- top_pfams %>%
  left_join(df %>% select(Pfam1, DE1) %>% distinct(), by = "Pfam1") %>%
  mutate(label = paste(Function, Pfam1, DE1, sep = " | "))

pfam_label_map <- setNames(label_df$label, label_df$Pfam1)
names(pfam_norm)[-1] <- pfam_label_map[names(pfam_norm)[-1]]

# === Step 8: Melt into long format ===
heatmap_long <- pfam_norm %>%
  pivot_longer(-Order, names_to = "Pfam_Label", values_to = "Percent")

# Convert to character to safely modify
heatmap_long$Pfam_Label <- as.character(heatmap_long$Pfam_Label)

# Shorten specific labels
heatmap_long$Pfam_Label[heatmap_long$Pfam_Label == "Signaling | PF04283 | Chemotaxis signal transduction system protein F from archaea"] <- 
  "Signaling | PF04283 | Archaeal signal transduction protein CheF"

heatmap_long$Pfam_Label[heatmap_long$Pfam_Label == "Packaging | PF07141 | Putative bacteriophage terminase small subunit"] <- 
  "Packaging | PF07141 | Bacteriophage terminase small subunit"

heatmap_long$Pfam_Label[heatmap_long$Pfam_Label == "Transport | PF03364 | Polyketide cyclase / dehydrase and lipid transport"] <- 
  "Transport | PF03364 | Polyketide cyclase/dehydrase and lipid transport"

# Reorder factor levels
heatmap_long$Pfam_Label <- factor(
  heatmap_long$Pfam_Label,
  levels = rev(sort(unique(heatmap_long$Pfam_Label)))
)

# Define breaks and labels for heatmap bins
value_breaks <- c(-Inf, 0, 0.01, 0.05, 0.15, 0.50, 0.85, 0.99, Inf)
labels <- c("0", "<1%", "1–5%", "5–15%", "15–50%", "50–85%", "85–99%", "Core (99%+)")


# Bin values into categories
heatmap_long$value_category <- cut(
  heatmap_long$Percent, 
  breaks = value_breaks,
  labels = labels,
  include.lowest = TRUE
)

# Create color palette
colors <- brewer.pal(n = 9, name = "Blues")[2:9]
colors[1] <- "white"

# Final heatmap plot
p2 <- ggplot(heatmap_long, aes(x = Order, y = Pfam_Label, fill = value_category)) +
  geom_tile(color = "grey") +
  scale_fill_manual(values = colors, name = "% Genomes\nwith Pfam") +
  theme_minimal() +
  labs(
    title = "Most Conserved Pfam by Function Across Orders",
    x = "Order",
    y = "Function | Pfam | Definition"
  ) +
  theme(
    axis.text.x = element_text(size=12, angle = 90, hjust = 1, vjust = 0.4),
    axis.text.y = element_text(size = 12),
    axis.ticks.x = element_line(),
    plot.title = element_text(size = 18, hjust = 1),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12)
  )

p2


ggsave("../plots/heatmap_most_conserved_structures.tiff", plot = p, width = 8, height = 8, dpi = 600)

ggsave("../plots/heatmap_most_conserved_pfams.tiff", plot = p2, width = 8, height = 8, dpi = 600)


# Load library
library(patchwork)

# Combine vertically
combined_plot <- p / p2  # "/" stacks vertically
combined_plot

# Save to file (half A4 landscape width)
ggsave("../plots/combined_heatmaps_vertical.tiff", 
       plot = combined_plot, 
       width = 9.5, 
       height = 11.69,
       units = "in", 
       dpi = 600)



