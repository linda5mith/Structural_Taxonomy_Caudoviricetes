library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(dplyr)

# Read the data
setwd('C:/Users/crtuser/Documents/Structural_Taxonomy/data')

pairwise_df <- read.csv('4083_%_shared_proteins.csv', row.names = 1)
RED_taxonomy <- read.csv('RED_taxonomy.csv',encoding = "UTF-8")

# Clean taxonomy to make sure no NAs
RED_taxonomy_clean <- RED_taxonomy[!(is.na(RED_taxonomy$Family) | RED_taxonomy$Family == ""), ]

# Count occurrences of each Family and filter out those with ≤ 8 members
RED_taxonomy_filtered <- RED_taxonomy_clean %>%
  group_by(Family) %>%
  filter(n() >= 100) %>%
  ungroup()

dim(RED_taxonomy_filtered)

# Extract the valid RefSeq accessions (Leaves)
valid_leaves <- RED_taxonomy_filtered$Leaves

# Now, align the rownames and colnames
pairwise_subset_filtered <- pairwise_df[valid_leaves, valid_leaves]

print(dim(pairwise_subset_filtered))

#pairwise_subset <- pairwise_subset_filtered[1:2000, 1:2000]
pairwise_matrix <- as.matrix(pairwise_subset_filtered)

# Get family information for the subset based on Virus.Refseq.Accession
family_info <- setNames(RED_taxonomy_filtered$Family, RED_taxonomy_filtered$Leaves)

# Extract family names for rows and columns of the pairwise subset
row_families <- family_info[rownames(pairwise_matrix)]
col_families <- family_info[colnames(pairwise_matrix)]

print(unique(row_families))
print(unique(col_families))

# Check are any row names not found in family_info
missing_rows <- rownames(pairwise_matrix)[!rownames(pairwise_matrix) %in% names(family_info)]
print(missing_rows)  # Should be 0

################### Defining colours ################################

# Get the unique families from the rows and columns of pairwise_matrix
unique_row_families <- unique(family_info[rownames(pairwise_matrix)])
unique_col_families <- unique(family_info[colnames(pairwise_matrix)])

# Combine both row and column families, and get the unique ones
unique_families <- unique(c(unique_row_families, unique_col_families))
print(paste('Unique families..', length(unique(unique_families))))

# Generate 138 unique colors
# Set seed for reproducibility
set.seed(10)
library(wesanderson)

# Define the family colors (create a palette)
palette1 <- brewer.pal(8, "Set3")
#palette2 <- brewer.pal(8, "Set1")
palette3 <- brewer.pal(8, "Dark2")
palette4 <- viridis(12)
palette5 <- viridis(12, option = "D")
palette6 <- viridis(12, option = "C")
palette7 <- brewer.pal(9, "Spectral")
palette8 <- brewer.pal(9, "RdYlBu")
palette9 <- wes_palette("Zissou1", 5, type = "discrete")
palette10 <- wes_palette("GrandBudapest1", 4, type = "discrete")
#palette11 <- magma(12)
palette12 <- inferno(12)
palette13 <- cividis(12)
palette14 <- wes_palette("Moonrise3", 5, type = "discrete")
palette15 <- wes_palette("GrandBudapest2", 4, type = "discrete")
palette16 <- wes_palette("Rushmore1", 5, type = "discrete")
palette17 <- wes_palette("Royal2", 5, type = "discrete")
palette18 <- wes_palette("Darjeeling2", 5, type = "discrete")
palette19 <- wes_palette("FantasticFox1", 5, type = "discrete")
palette20 <- wes_palette("Zissou1", 5, type = "discrete")
palette21 <- wes_palette("Chevalier1", 4, type = "discrete")

combined_palette <- c(palette1, palette3, palette4, palette5, palette6, palette7,
                      palette8, palette9, palette10,palette12,palette13,palette14,palette15,palette16,
                      palette17,palette18,palette19, palette20, palette21)

print(paste('Combined palette length..', length(unique(combined_palette))))

# # Add more colors from the 'Paired' palette
# additional_colors <- brewer.pal(12, "Paired")
# combined_palette2 <- c(combined_palette, additional_colors)
sum(is.na(RED_taxonomy_filtered$Family))  # should be 0

# Calculate the number of unique families (ensure no duplicates)
unique_families <- unique(RED_taxonomy_filtered$Family)

# Make sure you have enough colors for unique families
additional_needed = length(unique_families) - length(combined_palette)
print(additional_needed)  # ✅ this prints the number of additional colors needed

# Generate additional unique colors if necessary
if (additional_needed > 0) {
  print(paste(additional_needed, "Additional colors needed.."))
  extra_colors <- colorspace::rainbow_hcl(additional_needed, start = 120, end = 240)  # Only blues & greens
  combined_palette <- unique(c(combined_palette, extra_colors))
}

# Ensure no duplicates in the final palette
final_palette <- unique(combined_palette)
print(paste('Final palette length..', length(unique(final_palette))))

# Ensure we pad it to the correct length manually
while (length(final_palette) < length(unique_families)) {
  final_palette <- c(final_palette, generate_random_color())
  final_palette <- unique(final_palette)
}


# (Optional) Display the palette
barplot(rep(1, length(final_palette)), col=final_palette, border=NA, main="Unique Colors for Families")


#final_palette <- sample(final_palette)
# Assign each family a unique color
family_colors <- setNames(final_palette[1:length(unique_families)], unique_families)

# Function to detect dark or black colors
is_dark <- function(color) {
  rgb_vals <- col2rgb(color) / 255  # Convert hex to RGB values between 0 and 1
  luminance <- 0.2126 * rgb_vals[1] + 0.7152 * rgb_vals[2] + 0.0722 * rgb_vals[3]  # Calculate luminance
  return(luminance < 0.2)  # Dark colors have low luminance
}

# Generate a random color function
generate_random_color <- function() {
  rgb(runif(1), runif(1), runif(1))  # Random RGB color
}

# Replace dark colors with random colors
family_colors_fixed <- family_colors
family_colors_fixed[which(sapply(family_colors, is_dark))] <- sapply(1:sum(sapply(family_colors, is_dark)), function(x) generate_random_color())

#family_colors_fixed <- sample(family_colors_fixed)

# Assign each family a unique color
family_colors <- setNames(family_colors_fixed[1:length(unique_families)], unique_families)
barplot(rep(1, length(family_colors_fixed)), col=family_colors_fixed, border=NA, main="Unique Colors for Families")

family_colors_fixed <- sample(family_colors_fixed)

family_colors <- family_colors_fixed
barplot(rep(1, length(family_colors)), col=family_colors, border=NA, main="Unique Colors for Families")


#################### Initial heatmap #############################

# Split families for row and column annotations
row_family_split <- factor(row_families)
col_family_split <- factor(col_families)

# Heatmap ranges colours
colors <- brewer.pal(n = 9, name = "YlGnBu")

# Manually add "white" as the first color, then interpolate to reach 12 colors
color_palette <- colorRampPalette(c("white", colors))(12)

memory.limit(size = 56000)  # Increase memory limit to 56GB

initial_heatmap <- Heatmap(
  pairwise_matrix,
  col = colorRamp2(
    c(0, 0.00001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
    color_palette
  ),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_dend_gp = gpar(lwd = 0.3),
  column_dend_gp = gpar(lwd = 0.3),
  column_title = "% Shared Structural Proteins by Family",
  column_title_gp = gpar(fontsize = 25),
  name = '% proteins\n shared',
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_split = row_family_split,
  column_split = col_family_split,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 4),
  use_raster = TRUE,
  raster_device = "png",  # Force PNG rasterization
  raster_quality = 3,  # Reduce quality slightly to avoid overload
  row_gap = unit(0.7, "mm"),
  column_gap = unit(0.7, "mm")
)

#pdf("perc_shared_heatmap_INITIAL_03_03_25.pdf", width = 18, height = 18)
ht <- draw(initial_heatmap)
#dev.off()

# # Save the heatmap
# png("perc_shared_heatmap_INITIAL_01_03_25.png", width = 16, height = 16, units = "in", res = 700)
# ht <- draw(initial_heatmap)
# dev.off()

# Get row order without extra indexing
ordered_row_indices <- unlist(row_order(ht), use.names = FALSE)
ordered_row_families <- row_families[ordered_row_indices]

# Get column order similarly if needed
ordered_col_indices <- unlist(column_order(ht), use.names = FALSE)
ordered_col_families <- col_families[ordered_col_indices]

# Verify the ordered families
print(ordered_row_families)
print(ordered_col_families)

# Match order of columns and rows to family-assigned colours
col_families_unique <- unique(ordered_col_families)
row_families_unique <- unique(ordered_row_families)

family_colors_fixed <- sample(family_colors_fixed)

# Define a new palette based on the order of col_families_unique
col_colours <- setNames(family_colors_fixed[1:length(col_families_unique)], col_families_unique)
row_colours <- setNames(family_colors_fixed[1:length(row_families_unique)], col_families_unique)

# Check that the colors have been assigned correctly
print(col_colours)

# If you want to use this for a plot, for example, a bar plot:
barplot(rep(1, length(col_colours)), col = col_colours, border = NA, main = "Unique Colors for Column Families")

# Convert to data frame
colours_df <- data.frame(
  Family = names(col_colours),
  Colour = unname(col_colours),
  stringsAsFactors = FALSE
)


col_colours <- setNames(colours_df$Colour, colours_df$Family)

col_colours <- family_colors[match(col_families_unique, names(family_colors))]
row_colours <- family_colors[match(row_families_unique, names(family_colors))]


newly_defined_color <- "blue"   # Color for families containing "Family"
ictv_vmr40_color <- "black"     # Color for standard ICTV families
newly_defined_type <- "dashed"  # Dashed border for new families
ictv_vmr40_type <- "solid"      # Solid border for ICTV families

# Create block border color mapping for columns
matched_indices <- match(col_families_unique, RED_taxonomy_filtered$Family)
matched_definitions <- RED_taxonomy_filtered$Family_Source[matched_indices]

# Assign colors based on Family.Definition
family_border_col_colors <- ifelse(matched_definitions == "ICTV VMR40", ictv_vmr40_color, newly_defined_color)
family_border_col_types <- ifelse(matched_definitions == "ICTV VMR40", ictv_vmr40_type, newly_defined_type)

# Create block border color mapping for rows
family_border_row_colors <- ifelse(matched_definitions == "ICTV VMR40", ictv_vmr40_color, newly_defined_color)
family_border_row_types <- ifelse(matched_definitions == "ICTV VMR40", ictv_vmr40_type, newly_defined_type)

# Print results
print(family_border_col_colors)
print(family_border_col_types)
print(family_border_row_colors)
print(family_border_row_types)
# # Column Annotation
top_family_annotation <- HeatmapAnnotation(
  Family = anno_block(
    gp = gpar(
      fill = col_colours,
      col = family_border_col_colors,
      lwd = 0.2,
      lty = family_border_col_types
    ),
    labels = unique(ordered_col_families),
    labels_gp = gpar(col = "black", fontsize = 5),
    labels_rot = 90
  )
)
# 
# # Row Annotation
side_family_annotation <- rowAnnotation(
  Family = anno_block(
    gp = gpar(
      fill = row_colours,
      col = family_border_row_colors,
      lwd = 0.2,
      lty = family_border_row_types
    ),
    labels = unique(ordered_row_families),
    labels_gp = gpar(col = "black", fontsize = 5),
    labels_rot = 0  # Rotating row labels for better readability
  )
)

newly_defined_color <- "#5C7285"  
ictv_vmr40_color2 <- "#C5D3E8"  

# Fix family definition type assignment (ensure correct factor levels)
family_bottom_block_colours <- ifelse(matched_definitions == "ICTV VMR40", ictv_vmr40_color2, newly_defined_color)

# Define explicit color mapping
family_colors_def <- c("Newly Defined" = newly_defined_color, "ICTV VMR40" = ictv_vmr40_color2)

# Create bottom family annotation with correct color mapping
bottom_family_annotation <- HeatmapAnnotation(
  Family_Definition = anno_block(
    gp = gpar(
      fill = family_bottom_block_colours,  # Ensure correct color mapping
      col = "black",
      lwd = 0.1
    ),
    labels_gp = gpar(col = "black", fontsize = 6),
    which = "column"  
  ),
  annotation_name_gp = gpar(col = NA)  # Hide annotation name
)


ht <- Heatmap(
  pairwise_matrix,
  col = colorRamp2(
    c(0, 0.00001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
    color_palette),
  cluster_rows = TRUE,       
  cluster_columns = TRUE,    
  row_dend_gp = gpar(lwd = 0.2),
  column_dend_gp = gpar(lwd = 0.2),
  column_title = "% Shared Structural Proteins by Family",
  column_title_gp = gpar(fontsize = 30),
  name = '% proteins\n shared',
  show_heatmap_legend = FALSE,  
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_split = row_family_split,
  column_split = col_family_split,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 4),
  top_annotation = top_family_annotation,
  left_annotation = side_family_annotation,
  bottom_annotation = bottom_family_annotation,
  use_raster = TRUE,
  raster_device = "png",  # Force PNG rasterization
  raster_quality = 2,  # Reduce quality slightly to avoid overload
)

# Manually extract the color mapping legend to avoid NULL errors
ht_legend <- color_mapping_legend(ht@matrix_color_mapping, title = "% Proteins Shared")

# Define the new legend for family definition
legend_family_definition <- Legend(
  labels = names(family_colors_def),
  legend_gp = gpar(fill = family_colors_def),
  title = "Family Definition"
)

# Stack the legends together
combined_legends <- packLegend(ht_legend, legend_family_definition, direction = "vertical")
draw(combined_legends, x = unit(1, "npc") - unit(4, "mm"), just = "right")

# Save the heatmap
pdf("../plots/perc_shared_proteins_by_family_heatmap.pdf", width = 18, height = 18)
draw(ht, annotation_legend_list = list(combined_legends), padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()


# Open the PNG device
png("../plots/perc_shared_proteins_by_family_heatmap.png", width = 16, height = 16, units = "in", res = 700)
draw(ht, annotation_legend_list = list(combined_legends), padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()

