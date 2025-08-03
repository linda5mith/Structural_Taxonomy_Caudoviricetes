library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(dplyr)



set.seed(123)
mat2 = matrix(rnorm(50*50), nrow = 50)
split = rep(1:5, each = 10)

ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  foo = anno_block(gp = gpar(fill = 2:6), labels = LETTERS[1:5])
)
Heatmap(mat2, name = "mat2", column_split = split, top_annotation = ha, 
        column_title = NULL)

library(GetoptLong)  # for the function qq()
group_block_anno = function(group, empty_anno, gp = gpar(), 
                            label = NULL, label_gp = gpar()) {
  
  seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
  loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
  loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  
  seekViewport("global")
  grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
            just = c("left", "bottom"), gp = gp)
  if(!is.null(label)) {
    grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = label_gp)
  }
}

group_block_anno(1:3, "empty", gp = gpar(fill = "red"), label = "group 1")
group_block_anno(4:5, "empty", gp = gpar(fill = "blue"), label = "group 2")


################################################################


###############################################################




# Read the data
setwd('C:/Users/crtuser/Documents/CrAss PHD/Crassvirales_ICTV_proposal_2024/Structural_Taxonomy/data')

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

# ========= Generate Order Colors ==========
most_common_families <- RED_taxonomy_filtered %>%
  filter(!is.na(Order)) %>%
  group_by(Order, Family) %>%
  tally() %>%
  arrange(Order, desc(n)) %>%
  slice_head(n = 1) %>%
  ungroup()
order_colors <- setNames(family_colors[most_common_families$Family], most_common_families$Order)

# ========= Heatmap Annotations ==========
column_split <- factor(col_families)
names(column_split) <- colnames(pairwise_matrix)
row_split <- factor(row_families)

# Family-level annotations
col_families_unique <- unique(col_families)
row_families_unique <- unique(row_families)
col_colours <- family_colors[col_families_unique]
row_colours <- family_colors[row_families_unique]

matched_indices <- match(col_families_unique, RED_taxonomy_filtered$Family)
matched_definitions <- RED_taxonomy_filtered$Family_Source[matched_indices]
family_border_col_colors <- ifelse(matched_definitions == "ICTV VMR40", "black", "blue")
family_border_col_types <- ifelse(matched_definitions == "ICTV VMR40", "solid", "dashed")


family_top <- HeatmapAnnotation(
  empty_order_col = anno_empty(border = FALSE, height = unit(4, "mm")),  # NEW
  # empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  Family = anno_block(
    gp = gpar(fill = col_colours, col = family_border_col_colors, lwd = 0.2, lty = family_border_col_types),
    labels = col_families_unique,
    labels_gp = gpar(col = "black", fontsize = 5),
    labels_rot = 90
  ),
  annotation_name_side = "left",
  which = "column"
)

# # Row Annotation
side_family_annotation <- rowAnnotation(
  empty_order_row = anno_empty(border = FALSE, height = unit(8, "mm")),
  Family = anno_block(
    gp = gpar(
      fill = row_colours,
      col = family_border_col_colors,
      lwd = 0.2,
      lty = family_border_col_types
    ),
    labels = row_families_unique,
    labels_gp = gpar(col = "black", fontsize = 5),
    labels_rot = 0  # Rotating row labels for better readability
  )
)


# ========= Final Heatmap ==========
color_palette <- colorRampPalette(c("white", brewer.pal(9, "YlGnBu")))(12)
ht <- Heatmap(
  pairwise_matrix,
  col = colorRamp2(c(0, 0.00001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), color_palette),
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
  row_split = row_split, 
  column_split = column_split,
  # remove spacing between blocks
  row_gap = unit(0, "mm"),
  column_gap = unit(0, "mm"),
  row_title_rot = 0, 
  row_title_gp = gpar(fontsize = 4),
  top_annotation = family_top,
  left_annotation = side_family_annotation,
  use_raster = TRUE, 
  raster_device = "png", raster_quality = 2
)


# Draw heatmap and add side annotation
#png("../plots/perc_shared_proteins_with_order_blocks.png", width = 16, height = 16, units = "in", res = 700)
ht_drawn <- draw(ht, 
                 padding = unit(c(5, 5, 5, 5), "mm"))

# Get row order without extra indexing
ordered_row_indices <- unlist(row_order(ht), use.names = FALSE)
ordered_row_families <- row_families[ordered_row_indices]

# Get column order similarly if needed
ordered_col_indices <- unlist(column_order(ht), use.names = FALSE)
ordered_col_families <- col_families[ordered_col_indices]

# Verify the ordered families
#print(ordered_row_families)
#print(ordered_col_families)


group_block_anno <- function(group, empty_anno, label,
                             gp = gpar(), label_gp = gpar(fontsize = 6),
                             horizontal = TRUE) {
  seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
  loc1 <- deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  
  seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
  loc2 <- deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  
  seekViewport("global")
  
  if (horizontal) {
    # Column direction
    grid.rect(loc1$x, loc1$y,
              width = loc2$x - loc1$x,
              height = unit(5, "mm"),
              just = c("left", "bottom"),
              gp = gp)
    
    grid.text(label,
              x = (loc1$x + loc2$x) / 2,
              y = loc1$y + unit(3.2, "mm"),
              just = "bottom",
              gp = label_gp)
    
  } else {
    # Row direction
    grid.rect(loc1$x, loc1$y,
              width = unit(3, "mm"),
              height = loc2$y - loc1$y,
              just = c("left", "bottom"),
              gp = gp)
    
    grid.text(label,
              x = loc1$x + unit(3.2, "mm"),
              y = (loc1$y + loc2$y) / 2,
              just = "left",
              gp = label_gp)
  }
}



#  Row Order Blocks by Family ####
rle_rows <- rle(ordered_row_families)
start_idx <- cumsum(c(1, head(rle_rows$lengths, -1)))
end_idx <- cumsum(rle_rows$lengths)
row_block_labels <- rle_rows$values  # These are Family names

row_orders_per_block <- family_to_order[row_block_labels]

unique_row_orders <- unique(row_orders_per_block)
row_order_colors <- structure(
  brewer.pal(length(unique_row_orders), "Set3")[seq_along(unique_row_orders)],
  names = unique_row_orders
)

#draw one order label per blocks...
rle_order_rows <- rle(ordered_row_orders)
row_start <- cumsum(c(1, head(rle_order_rows$lengths, -1)))
row_end <- cumsum(rle_order_rows$lengths)


for (i in seq_along(row_orders_per_block)) {
  this_order <- row_orders_per_block[i]
  
  if (!is.na(this_order)) {
    group_block_anno(
      group = i,
      empty_anno = "empty_order_row",
      label = this_order,
      gp = gpar(fill = row_order_colors[this_order], col = NA),
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
      gp = gpar(fill = row_order_colors[this_order], col = NA),
      label_gp = gpar(fontsize = 0, fontface = "bold"),
      horizontal = TRUE
    )
  }
}



col_orders_by_slice <- column_order(ht)
col_families_by_slice <- lapply(col_orders_by_slice, function(idxs) col_families[ordered_col_indices[idxs]])
col_orders_by_slice <- lapply(col_families_by_slice, function(fams) family_to_order[fams])


for (s in seq_along(col_orders_by_slice)) {
  order_vec <- col_orders_by_slice[[s]]
  rle_slice <- rle(order_vec)
  slice_start <- cumsum(c(1, head(rle_slice$lengths, -1)))
  slice_end <- cumsum(rle_slice$lengths)
  
  for (i in seq_along(rle_slice$values)) {
    this_order <- rle_slice$values[i]
    
    group_block_anno(
      group = slice_start[i]:slice_end[i],
      empty_anno = "empty_order_col",
      label = this_order,
      gp = gpar(fill = row_order_colors[this_order], col = NA),
      label_gp = gpar(fontsize = 6, fontface = "bold"),
      horizontal = TRUE
    )
  }
}




##################################


group_block_anno <- function(group, empty_anno, label,
                             gp = gpar(), label_gp = gpar(fontsize = 6),
                             horizontal = TRUE,
                             slice_id = 1) {
  seekViewport(qq("annotation_@{empty_anno}_@{slice_id}"))
  
  loc1 <- deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  loc2 <- deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  
  seekViewport("global")
  
  if (horizontal) {
    grid.rect(loc1$x, loc1$y,
              width = loc2$x - loc1$x,
              height = unit(5, "mm"),
              just = c("left", "bottom"),
              gp = gp)
    
    grid.text(label,
              x = (loc1$x + loc2$x) / 2,
              y = loc1$y + unit(2.5, "mm"),
              just = "center",
              gp = label_gp)
  } else {
    grid.rect(loc1$x, loc1$y,
              width = unit(3, "mm"),
              height = loc2$y - loc1$y,
              just = c("left", "bottom"),
              gp = gp)
    
    grid.text(label,
              x = loc1$x + unit(2.5, "mm"),
              y = (loc1$y + loc2$y) / 2,
              just = "center",
              gp = label_gp)
  }
}


for (s in seq_along(col_orders_by_slice)) {
  order_vec <- col_orders_by_slice[[s]]
  rle_slice <- rle(order_vec)
  slice_start <- cumsum(c(1, head(rle_slice$lengths, -1)))
  slice_end <- cumsum(rle_slice$lengths)
  
  for (i in seq_along(rle_slice$values)) {
    this_order <- rle_slice$values[i]
    
    group_block_anno(
      group = slice_start[i]:slice_end[i],  # relative within slice
      empty_anno = "empty_order_col",
      label = this_order,
      gp = gpar(fill = order_colors[this_order], col = NA),
      label_gp = gpar(fontsize = 6, fontface = "bold"),
      horizontal = TRUE
    )
    
    # Use correct viewport name
    seekViewport(qq("annotation_empty_order_col_@{s}"))
    # then use that context to calculate loc1, loc2
  }
}







###################################






library(GetoptLong)  # for `qq()` used in viewport names

# Map each family to its Order (adjust column names as needed)
family_to_order <- setNames(RED_taxonomy_filtered$Order, RED_taxonomy_filtered$Family)

# Define order groups based on the order of `ordered_col_families`
ordered_col_orders <- family_to_order[ordered_col_families]
ordered_row_orders <- family_to_order[ordered_row_families]

# Assign a numeric ID to each Order block
order_ids <- as.numeric(factor(ordered_col_orders))
order_labels <- levels(factor(ordered_col_orders))


# group order blocks function
group_block_anno <- function(group, empty_anno, label, gp = gpar(), label_gp = gpar(fontsize = 6)) {
  seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
  loc1 <- deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  
  seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
  loc2 <- deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  
  seekViewport("global")
  grid.rect(loc1$x, loc1$y,
            width = loc2$x - loc1$x,
            height = unit(3, "mm"),
            just = c("left", "bottom"),
            gp = gp)
  grid.text(label,
            x = (loc1$x + loc2$x) / 2,
            y = loc1$y + unit(3.2, "mm"),
            just = "bottom",
            gp = label_gp)
}



# Get column order by slice
col_orders_by_slice <- column_order(ht)
col_families_by_slice <- lapply(col_orders_by_slice, function(idxs) col_families[ordered_col_indices[idxs]])
col_orders_by_slice <- lapply(col_families_by_slice, function(fams) family_to_order[fams])

# Define colors for orders
order_colors <- structure(
  brewer.pal(length(unique(unlist(col_orders_by_slice))), "Set3")[seq_along(unique(unlist(col_orders_by_slice)))],
  names = unique(unlist(col_orders_by_slice))
)

# Iterate through each slice (annotation_empty_order_1 to _N)
for (i in seq_along(col_orders_by_slice)) {
  current_orders <- col_orders_by_slice[[i]]
  current_order <- unique(current_orders[!is.na(current_orders)])
  
  if (length(current_order) == 1) {
    group_block_anno(
      group = i,
      empty_anno = "empty_order_col",
      label = current_order,
      gp = gpar(fill = order_colors[current_order], col = NA),
      label_gp = gpar(fontsize = 4, fontface = "bold")
    )
  } else if (length(current_order) > 1) {
    # Optional: if slice contains mixed Orders (shouldn't happen), collapse label
    group_block_anno(
      group = i,
      empty_anno = "empty_order_col",
      label = paste(current_order, collapse = ", "),
      gp = gpar(fill = "grey", col = NA),
      label_gp = gpar(fontsize = 6, fontface = "italic")
    )
  }
}


### Draw side annotation order blocks
row_orders_by_slice <- row_order(ht)
row_families_by_slice <- lapply(row_orders_by_slice, function(idxs) row_families[ordered_row_indices[idxs]])
row_orders_by_slice <- lapply(row_families_by_slice, function(fams) family_to_order[fams])

# Use same `order_colors` from above

for (i in seq_along(row_orders_by_slice)) {
  current_orders <- row_orders_by_slice[[i]]
  current_order <- unique(current_orders[!is.na(current_orders)])
  
  if (length(current_order) == 1) {
    group_block_anno(
      group = i,
      empty_anno = "empty_order_row",
      label = current_order,
      gp = gpar(fill = order_colors[current_order], col = NA),
      label_gp = gpar(fontsize = 6, fontface = "bold"),
      horizontal = FALSE
    )
  } else if (length(current_order) > 1) {
    group_block_anno(
      group = i,
      empty_anno = "empty_order_row",
      label = paste(current_order, collapse = ", "),
      gp = gpar(fill = "grey", col = NA),
      label_gp = gpar(fontsize = 5, fontface = "italic"),
      horizontal = FALSE
    )
  }
}




###########################################################
# ========= Add Order Block Layer ==========
col_ordered <- colnames(pairwise_matrix)[unlist(column_order(ht_drawn))]
ordered_col_orders <- order_info[col_ordered]
rle_orders <- rle(ordered_col_orders)
start_idx <- cumsum(c(1, head(rle_orders$lengths, -1)))
end_idx <- cumsum(rle_orders$lengths)



for (i in seq_along(rle_orders$values)) {
  ord <- rle_orders$values[i]
  if (!is.na(ord)) {
    group_block_anno(start_idx[i]:end_idx[i], "empty", ord, order_colors[ord])
  }
}

#dev.off()


