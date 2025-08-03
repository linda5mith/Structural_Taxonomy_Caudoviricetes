library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(rlang)

# Read the data
setwd('C:/Users/crtuser/Documents/CrAss PHD/Crassvirales_ICTV_proposal_2024/Structural_Taxonomy/data')

pairwise_df <- read.csv('4083_%_shared_proteins.csv', row.names = 1)
RED_taxonomy <- read.csv('RED_taxonomy.csv',encoding = "UTF-8")

# Clean taxonomy to make sure no NAs
RED_taxonomy_clean <- RED_taxonomy[!(is.na(RED_taxonomy$Family) | RED_taxonomy$Family == ""), ]

# Count occurrences of each Family and filter out those with ≤ 8 members
RED_taxonomy_filtered <- RED_taxonomy_clean %>%
  group_by(Family) %>%
  filter(n() >= 10) %>%
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
combined_palette <- c(
  brewer.pal(8, "Set3"),
  brewer.pal(8, "Dark2"),
  viridis(12),
  viridis(12, option = "D"),
  viridis(12, option = "C"),
  brewer.pal(9, "Spectral"),
  brewer.pal(9, "RdYlBu"),
  wes_palette("Zissou1", 5, type = "discrete"),
  wes_palette("GrandBudapest1", 4, type = "discrete"),
  inferno(12),
  cividis(12),
  wes_palette("Moonrise3", 5, type = "discrete"),
  wes_palette("GrandBudapest2", 4, type = "discrete"),
  wes_palette("Rushmore1", 5, type = "discrete"),
  wes_palette("Royal2", 5, type = "discrete"),
  wes_palette("Darjeeling2", 5, type = "discrete"),
  wes_palette("FantasticFox1", 5, type = "discrete"),
  wes_palette("Zissou1", 5, type = "discrete"),
  wes_palette("Chevalier1", 4, type = "discrete")
)


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


colours_df <- read.csv('../../STRUCTURAL_TAXONOMY_04_03_25/data/family_heatmap_colours2.csv')

# 2. Make a named color vector: names are families, values are hex colors
family_colors <- setNames(colours_df$Colour, colours_df$Family)

# 3. Extract colors for column and row families
col_colours <- family_colors[match(col_families_unique, names(family_colors))]
row_colours <- family_colors[match(row_families_unique, names(family_colors))]


# col_colours <- family_colors[match(col_families_unique, names(family_colors))]
# row_colours <- family_colors[match(row_families_unique, names(family_colors))]


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


#######################################################################################


## === Define Family → Order Mapping ===
family_to_order <- setNames(RED_taxonomy_filtered$Order, RED_taxonomy_filtered$Family)

## === Define Dominant Family per Order ===
family_counts <- RED_taxonomy_filtered %>%
  count(Order, Family, sort = TRUE)

dominant_family_per_order <- family_counts %>%
  group_by(Order) %>%
  top_n(1, n) %>%
  ungroup() %>%
  select(Order, Family) %>%
  as.data.frame() %>%
  { setNames(.$Family, .$Order) }

## === Define Order Colors Based on Dominant Families ===
row_order_colors <- family_colors[dominant_family_per_order]
names(row_order_colors) <- names(dominant_family_per_order)

## === Define Annotations ===

# Top (Column) Family Annotation
top_family_annotation <- HeatmapAnnotation(
  empty_order_col = anno_empty(border = FALSE, height = unit(4, "mm")),
  Family = anno_block(
    height = unit(18, "mm"),  # <----- increase block height here
    gp = gpar(
      fill = col_colours,
      #col = family_border_col_colors,
      col = 'black',
      fontface = "bold",
      lwd = 0.2 #,
      #lty = family_border_col_types
    ),
    labels = unique(ordered_col_families),
    labels_gp = gpar(col = "black", 
                     fontsize = 5),
    labels_rot = 90
  ),
  annotation_label = c(" ", "Family"),  # label right of block
  annotation_name_side = "left",        # <<<<<< ADD THIS
  annotation_name_gp = gpar(fontsize = 6, fontface = "bold", col = 'black'),
  annotation_name_rot = 0               # horizontal label
)


# Left (Row) Family Annotation
side_family_annotation <- rowAnnotation(
  #empty_order_row = anno_empty(border = FALSE, width = unit(4, "mm")),
  Family = anno_block(
    width = unit(18, "mm"),  # <----- increase block width here
    gp = gpar(
      fill = row_colours,
      col = 'black',
      #col = family_border_row_colors,
      lwd = 0.2 #,
      #lty = family_border_row_types
    ),
    labels = unique(ordered_row_families),
    # text labels inside block
    labels_gp = gpar(col = "black", 
                     fontsize = 5),
    labels_rot = 0  # Rotating row labels for better readability
  )
)

# Bottom Family Definition Annotation
family_colors_def <- c(
  "Newly Defined" = "#5C7285",
  "ICTV VMR40" = "#C5D3E8"
)
family_bottom_block_colours <- ifelse(
  matched_definitions == "ICTV VMR40",
  family_colors_def["ICTV VMR40"],
  family_colors_def["Newly Defined"]
)

bottom_family_annotation <- HeatmapAnnotation(
  Family_Definition = anno_block(
    gp = gpar(fill = family_bottom_block_colours, col = "black", lwd = 0.1),
    labels_gp = gpar(col = "black", fontsize = 6),
    which = "column"
  ),
  annotation_name_gp = gpar(col = NA)
)


## === Define Heatmap ===
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
  # remove side annot dendro labels
  row_title = NULL,
  # remove spacing between blocks
  row_gap = unit(0, "mm"),
  column_gap = unit(0, "mm"),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 4),
  top_annotation = top_family_annotation,
  left_annotation = side_family_annotation,
  bottom_annotation = bottom_family_annotation,
  use_raster = TRUE,
  raster_device = "png",  # Force PNG rasterization
  raster_quality = 2,  # Reduce quality slightly to avoid overload
)

## === Legends ===
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



## === Order Block Drawing Function ===
# modified function to add order name to blocks...
## === Order Block Drawing Function ===
# modified function to add order name to blocks...
draw_order_blocks_from_families <- function(ordered_families,
                                            family_to_order,
                                            annotation_name = "Family",
                                            order_colors,
                                            axis = c("row", "column"),
                                            label_gp = gpar(fontsize = 6, fontface = "bold")) {
  axis <- match.arg(axis)
  
  ordered_orders <- family_to_order[ordered_families]
  rle_orders <- rle(ordered_orders)
  starts <- cumsum(c(0, head(rle_orders$lengths, -1))) + 1.5
  ends <- cumsum(rle_orders$lengths) + 0.5
  orders <- rle_orders$values
  
  for (i in seq_along(orders)) {
    this_order <- orders[i]
    if (!is.na(this_order)) {
      color <- if (this_order %in% names(order_colors)) order_colors[[this_order]] else "#CCCCCC"
      
      if (axis == "column") {
        decorate_annotation(annotation_name, {
          mid_x <- unit((starts[i] + ends[i] - 1) / 2, "native")
          grid.rect(
            x = mid_x,
            y = unit(1, "npc"),
            width = unit(ends[i] - starts[i] + 1, "native"),
            height = unit(5, "mm"),
            just = "bottom",
            gp = gpar(fill = color, col = "black", lwd = 0.2)
          )
          # grid.text(
          #   label = this_order,
          #   x = mid_x,
          #   y = unit(1, "npc") + unit(2.5, "mm"),  
          #   just = "center",                          
          #   gp = label_gp
          # )
        })
      } else {
        decorate_annotation(annotation_name, {
          mid_y <- unit((starts[i] + ends[i] - 1) / 2, "native")
          grid.rect(
            x = unit(0, "npc") - unit(3, "mm"),
            y = mid_y,
            width = unit(4.5, "mm"),
            height = unit(ends[i] - starts[i] + 1, "native"),
            just = "right",
            gp = gpar(fill = color, col = "black", lwd = 0.5)
          )
          grid.text(
            label = this_order,
            x = unit(0, "npc") - unit(5, "mm"),
            y = mid_y,
            just = "center",
            gp = label_gp
          )
        })
      }
    }
  }
}


## === Output to PDF ===
pdf("../data/plots/perc_shared_proteins_by_family_heatmap_order_blocks.pdf", width = 14, height = 14)
#png("../data/plots/perc_shared_proteins_by_family_heatmap_order_blocks.png", 
 #   width = 3000, height = 3000, res = 300)


ht_drawn <- draw(ht, annotation_legend_list = list(combined_legends), padding = unit(c(10, 10, 10, 10), "mm"))

# Draw order blocks after heatmap is drawn
draw_order_blocks_from_families(
  ordered_families = ordered_col_families,
  family_to_order = family_to_order,
  annotation_name = "Family",
  order_colors = row_order_colors,
  axis = "column",
  gpar(fontsize = 0, fontface = "bold", col = "white")
  #label_gp = gpar(fontsize = 6)
)

dev.off()

############################################################
############################################################

#### Ok so this code works :) just need to draw borders

for (i in seq_along(row_orders_per_block)) {
  this_order <- row_orders_per_block[i]
  
  if (!is.na(this_order)) {
    group_block_anno(
      group = i,
      empty_anno = "empty_order_row",
      label = this_order,
      gp = gpar(fill = row_order_colors[this_order], col = 'black'),
      label_gp = gpar(fontsize = 6, fontface = "bold"),
      horizontal = FALSE
    )
  }
}


for (i in seq_along(row_orders_per_block)) {
  this_order <- row_orders_per_block[i]
  
  if (!is.na(this_order)) {
    group_block_anno(
      group = i,
      empty_anno = "empty_order_col",
      label = this_order,
      gp = gpar(fill = row_order_colors[this_order], col = 'black'),
      label_gp = gpar(fontsize = 0, fontface = "bold"),
      horizontal = TRUE
    )
  }
}


#############################################

