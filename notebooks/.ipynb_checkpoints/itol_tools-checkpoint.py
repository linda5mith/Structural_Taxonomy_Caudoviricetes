import pandas as pd
import numpy as np
import textwrap
import random
import colorsys
from Bio import Phylo
import distinctipy
import colorcet as cc
import ast
import re
from matplotlib import cm
from matplotlib.colors import to_rgb, to_hex

# ------------------- COLOR UTILS -------------------

def generate_distinct_colors(n, seed=None):
    palette = list(cc.glasbey)
    if seed is not None:
        random.seed(seed)
        random.shuffle(palette)
    if n > len(palette):
        extra_needed = n - len(palette)
        print(f"🎨 Glasbey has {len(palette)} colors, generating {extra_needed} pastel fallback colors...")
        pastel_colors = [generate_pastel_color() for _ in range(extra_needed)]
        combo_palette = palette + pastel_colors
        random.shuffle(combo_palette)
        return combo_palette
    return palette

def generate_pastel_color():
    h, s, l = random.random(), 0.6, 0.7
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    return f'#{int(r*255):02X}{int(g*255):02X}{int(b*255):02X}'

def darken_color(hex_color, level=1, base_factor=0.9):
    """
    Darken a hex color progressively by multiplying RGB by base_factor^level.
    """
    rgb = to_rgb(hex_color)
    factor = base_factor ** level  # e.g., level 1 = 0.9, level 2 = 0.81, etc.
    dark_rgb = tuple(max(0, min(1, c * factor)) for c in rgb)
    return to_hex(dark_rgb)

# ------------------- FORMAT UTILS -------------------

def mrca_to_itol_format(s):
    """
    Converts MRCA like 'NC_000866_NC_048087_0.440762' to 'NC_000866|NC_048087'.
    Handles odd cases like:
    - 'BK063679_ON649698_0.473621' → 'BK063679|ON649698'
    - 'NC_041921_NC_041921' → 'NC_041921|NC_041921'
    """
    import re
    if not isinstance(s, str):
        print(f"❌ Invalid MRCA (not a string): {s}")
        return "INVALID_MRCA"

    matches = re.findall(r'[A-Z]{2}_[0-9]+|[A-Z]{2}[0-9]+', s)
    if len(matches) >= 2:
        return f"{matches[0]}|{matches[1]}"
    else:
        print(f"❌ Could not parse MRCA: {s}")
        return s
        
# ------------------- iTOL TREE COLOR -------------------

def create_itol_tree_colors_from_mrca(df, output_file, rank='Order', predefined_colors=None, seed=10):
    df = df.copy()
    df["Defined Rank Counts"] = df["Defined Rank Counts"].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

    label_list = []
    unique_labels = set()
    undefined_counter = 1
    for d in df["Defined Rank Counts"]:
        if not d:
            label = f"Undefined_{undefined_counter}"
            undefined_counter += 1
        else:
            sorted_counts = sorted(d.items(), key=lambda x: x[1], reverse=True)
            top_label = sorted_counts[0][0]
            if top_label.lower() == 'undefined':
                label = f"Undefined_{undefined_counter}"
                undefined_counter += 1
            else:
                label = top_label
        label_list.append(label)
        unique_labels.add(label)

    df["itol_clade"] = label_list

    if predefined_colors is None:
        predefined_colors = {}

    missing_labels = sorted(set(label_list) - set(predefined_colors.keys()))
    generated_colors = generate_distinct_colors(len(missing_labels), seed=seed)
    new_colors = dict(zip(missing_labels, generated_colors))
    full_color_mapping = {**predefined_colors, **new_colors}

    color_mapping = dict(zip(df["MRCA"], [full_color_mapping[label] for label in df["itol_clade"]]))
    label_mapping = dict(zip(df["MRCA"], df["itol_clade"]))

    header = f"""
    TREE_COLORS
    SEPARATOR SPACE
    SHOW_LEGEND 1
    LEGEND_TITLE {rank} Clade Colors
    DATA
    """

    output_lines = [textwrap.dedent(header.strip())]

    for (_, row) in df.iterrows():
        formatted = mrca_to_itol_format(row["MRCA"])
        clade_label = row['itol_clade']
        color = full_color_mapping[clade_label]
        output_lines.append(f"# {clade_label}")
        output_lines.append(f"{formatted} clade {color} normal 1")

    with open(output_file, 'w') as f:
        f.write("\n".join(output_lines))

    print(f"✅ TREE_COLORS file written to: {output_file}")
    return color_mapping, label_mapping

def create_itol_tree_colors_from_assigned_rank_v2(
    df,
    output_file,
    rank_col='New Order Name',
    legend_title='Clade Colors',
    color_seed=None,
    previous_annotation_file=None
):
    """
    Enhanced version: generates iTOL TREE_COLORS annotations, optionally reusing existing colors.
    """
    df = df.copy()
    df['itol_clade'] = df[rank_col].fillna('Undefined')
    unique_labels = df['itol_clade'].unique().tolist()

    if color_seed is not None:
        random.seed(color_seed)
        random.shuffle(unique_labels)

    # Reuse previous colors if provided
    previous_colors = {}
    used_colors = set()
    if previous_annotation_file:
        previous_colors = parse_previous_colors(previous_annotation_file)
        used_colors.update(previous_colors.values())

    label_to_color = {}
    new_color_index = 0
    generated_colors = generate_distinct_colors(len(unique_labels) + 10, seed=color_seed)

    for label in unique_labels:
        #base_label = label.split('_')[0]
        base_label = re.sub(r'_\d+$', '', label)
        if base_label in previous_colors:
            base_color = previous_colors[base_label]

            suffix_parts = label.split('_')
            try:
                level = int(suffix_parts[-1])
                if level == 1:
                    color = base_color
                    print(f"✅ '{label}' → base '{base_label}', using original color: {color}")
                else:
                    color = darken_color(base_color, level=level - 1)
                    print(f"🎨 '{label}' → darkened from '{base_label}', level={level - 1}: {base_color} → {color}")
            except ValueError:
                color = base_color
                print(f"✅ '{label}' → base '{base_label}', using original color (no numeric suffix): {color}")
        else:
            while generated_colors[new_color_index] in used_colors:
                new_color_index += 1
            color = generated_colors[new_color_index]
            new_color_index += 1
            print(f"🆕 New color for '{label}': {color}")

        label_to_color[label] = color
        used_colors.add(color)

    color_mapping = dict(zip(df['MRCA'], df['itol_clade'].map(label_to_color)))
    label_mapping = dict(zip(df['MRCA'], df['itol_clade']))

    header = f"""
    TREE_COLORS
    SEPARATOR SPACE
    SHOW_LEGEND 1
    LEGEND_TITLE {legend_title}
    DATA
    """
    output_lines = [textwrap.dedent(header.strip())]

    df = df.drop_duplicates(subset=['MRCA','Family'])

    for _, row in df.iterrows():
        mrca_fmt = mrca_to_itol_format(row['MRCA'])
        color = color_mapping.get(row['MRCA'], "#AAAAAA")
        output_lines.append(f"# {row['itol_clade']}")
        output_lines.append(f"{mrca_fmt} clade {color} normal 1")

    # Singleton branch coloring
    singleton_df = df[df['MRCA'].isna() | ~df['MRCA'].astype(str).str.contains(r'_', regex=True)]
    if not singleton_df.empty:
        for fam, subdf in singleton_df.groupby(rank_col):
            fam_color = label_to_color.get(fam, "#AAAAAA")
            output_lines.append(f"# Singleton family: {fam}")
            for leaf in subdf['Leaves']:
                if isinstance(leaf, str):
                    output_lines.append(f"{leaf} branch {fam_color} normal 1")

    # Encoding check
    for i, line in enumerate(output_lines):
        try:
            line.encode('cp1252')
        except UnicodeEncodeError:
            print(f"❌ Problem on line {i}: {repr(line)}")

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("\n".join(output_lines))

    print(f"✅ TREE_COLORS file written to: {output_file}")
    return color_mapping, label_mapping


# ------------------- MIDPOINT LEAF FINDER -------------------

def find_midpoint_leaf_from_mrca(tree, metadata_df, color_mapping, output_file, rank='Family', label='Clade labels', legend_title=None, label_mapping=None):
    if legend_title is None:
        legend_title = rank

    header = f"""
    DATASET_TEXT
    SEPARATOR COMMA

    #SHOW_LEGEND,1
    LEGEND_TITLE,{legend_title}

    DATASET_LABEL,{label}
    COLOR,#8174A0
    MARGIN,0
    SHOW_INTERNAL,0
    ALL_LABELS_ROTATION,0
    LABELS_BELOW,1
    VERTICAL_SHIFT,0
    STRAIGHT_LABELS,0
    ALIGN_TO_TREE,0
    SIZE_FACTOR,1
    EXTERNAL_LABEL_SHIFT,0

    DATA
    """

    output_lines = [textwrap.dedent(header.strip())]

    dedup_df = metadata_df.drop_duplicates(subset=['MRCA', rank])

    for _, row in dedup_df.iterrows():
        mrca_label = row.get('MRCA')

        clade_name = label_mapping.get(mrca_label, row.get(rank)) if label_mapping else row.get(rank)
        if pd.isna(mrca_label) or pd.isna(clade_name):
            continue

        clade_name = str(clade_name).strip().replace(',', '_')  # <-- this line is key

        color = color_mapping.get(mrca_label, "#000000")

        mrca_nodes = tree.search_nodes(name=mrca_label)
        if not mrca_nodes:
            print(f"⚠️ MRCA '{mrca_label}' not found in tree.")
            continue

        mrca_node = mrca_nodes[0]
        leaves = mrca_node.get_leaves()

        if not leaves:
            print(f"⚠️ No leaves under MRCA '{mrca_label}'.")
            continue

        total_dist = sum(leaf.dist for leaf in leaves)
        midpoint_leaf = None
        cumulative = 0
        for leaf in leaves:
            cumulative += leaf.dist
            if cumulative >= total_dist / 2:
                midpoint_leaf = leaf.name
                break

        if midpoint_leaf is None:
            print(f"⚠️ Could not determine midpoint leaf for MRCA '{mrca_label}'.")
            continue

        text_size = 7 if len(leaves) <= 3 else 17
        output_lines.append(f"{midpoint_leaf},{clade_name},-1,{color},bold,{text_size},0")

    with open(output_file, 'w') as f:
        f.write("\n".join(output_lines))

    print(f"✅ iTOL midpoint annotation file written to: {output_file}")

def find_midpoint_leaf_from_mrca_with_colorstrip(
    tree,
    metadata_df,
    colorstrip_file,
    output_file,
    rank='Order',
    label='Order labels',
    legend_title=None,
    label_mapping=None,
    mrca_col='MRCA',
    min_font=14,
    max_font=25,
    print_summary=True
):
    import textwrap
    import pandas as pd

    if legend_title is None:
        legend_title = rank

    # Load leaf → color mapping from colorstrip
    leaf_color_map = parse_colorstrip_file(colorstrip_file)

    header = f"""\
    DATASET_TEXT
    SEPARATOR COMMA

    LEGEND_TITLE,{legend_title}

    DATASET_LABEL,{label}
    COLOR,#8174A0
    MARGIN,0
    SHOW_INTERNAL,0
    ALL_LABELS_ROTATION,0
    LABELS_BELOW,1
    VERTICAL_SHIFT,0
    STRAIGHT_LABELS,0
    ALIGN_TO_TREE,0
    SIZE_FACTOR,1
    EXTERNAL_LABEL_SHIFT,0

    DATA
    """

    output_lines = [textwrap.dedent(header.strip())]
    seen_mrcas = set()
    labeled = []
    skipped = []
    not_found = []

    for _, row in metadata_df.iterrows():
        mrca_label = row.get(mrca_col)
        clade_name = label_mapping.get(mrca_label, row.get(rank)) if label_mapping else row.get(rank)

        if pd.isna(mrca_label) or pd.isna(clade_name):
            continue
        if mrca_label in seen_mrcas:
            continue
        seen_mrcas.add(mrca_label)

        clade_name = str(clade_name).strip().replace(',', '_')

        mrca_nodes = tree.search_nodes(name=mrca_label)
        if not mrca_nodes:
            not_found.append((mrca_label, clade_name))
            continue

        mrca_node = mrca_nodes[0]
        leaves = mrca_node.get_leaves()
        if not leaves:
            skipped.append((mrca_label, clade_name, "No leaves"))
            continue

        # Find midpoint leaf
        if len(leaves) == 1:
            midpoint_leaf = leaves[0].name
        else:
            total_dist = sum(leaf.dist for leaf in leaves)
            cumulative = 0
            midpoint_leaf = None
            for leaf in leaves:
                cumulative += leaf.dist
                if cumulative >= total_dist / 2:
                    midpoint_leaf = leaf.name
                    break
            if midpoint_leaf is None:
                midpoint_leaf = leaves[len(leaves) // 2].name  # fallback

        # Get color
        color = leaf_color_map.get(midpoint_leaf, "#000000")
        text_size = min_font if len(leaves) <= 5 else max_font

        if len(leaves) <= 5:
            print(f"🔍 Small clade '{clade_name}' ({len(leaves)} leaves) → font size {min_font} [midpoint: {midpoint_leaf}]")

        output_lines.append(f"{midpoint_leaf},{clade_name},-1,{color},bold,{text_size},0")
        labeled.append((mrca_label, clade_name, midpoint_leaf))

    # Write final file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("\n".join(output_lines))

    if print_summary:
        print(f"\n✅ iTOL midpoint label file written to: {output_file}")
        print(f"🟢 Labeled clades: {len(labeled)}")
        print(f"🔍 Small clades (<=5 leaves): {len([x for x in labeled if x[2] in leaf_color_map and leaf_color_map[x[2]] == '#000000'])}")
        print(f"⚠️ Skipped clades (no leaves): {len(skipped)}")
        print(f"❌ MRCAs not found in tree: {len(not_found)}")

        if not_found:
            print("\n❌ MRCAs not found in tree:")
            for mrca, name in not_found:
                print(f"  - {name} ({mrca})")

        if skipped:
            print("\n⚠️ Skipped clades (reason: no leaves):")
            for mrca, name, reason in skipped:
                print(f"  - {name} ({mrca}) → {reason}")


def create_itol_order_colorstrip_from_existing_colors(df, existing_colorstrip, output_file):
    # 1. Parse existing colorstrip (leaf → color)
    leaf_to_color = parse_colorstrip_file(existing_colorstrip)

    # 2. Get Order → most common color based on existing leaf colors
    df['Order'] = df['Order'].str.strip()
    df['Leaves'] = df['Leaves'].str.strip()
    
    order_color_map = {}
    for order, group in df.groupby("Order"):
        colors = group["Leaves"].map(leaf_to_color).dropna()
        if not colors.empty:
            most_common_color = colors.value_counts().idxmax()
            order_color_map[order] = most_common_color

    # 3. Assign new colors to orders not in original
    missing_orders = [o for o in df['Order'].unique() if o not in order_color_map]
    if missing_orders:
        import distinctipy
        new_colors = distinctipy.get_colors(len(missing_orders))
        hex_colors = [distinctipy.get_hex(c) for c in new_colors]
        for order, color in zip(missing_orders, hex_colors):
            order_color_map[order] = color
        print(f"🎨 Assigned new colors to {len(missing_orders)} orders")

    # 4. Build iTOL file
    header = """\
    DATASET_COLORSTRIP
    SEPARATOR COMMA
    DATASET_LABEL,Order_Annotations
    COLOR,#ff0000
    COLOR_BRANCHES,0
    BORDER_WIDTH,0.5
    BORDER_COLOR,#000000
    SHOW_STRIP_LABELS,0
    DATA
    """
    lines = [textwrap.dedent(header.strip())]
    for _, row in df.iterrows():
        leaf = row['Leaves']
        order = row['Order']
        color = order_color_map.get(order)
        if color:
            lines.append(f"{leaf},{color},")
        else:
            print(f"⚠️ Missing color for order '{order}'")

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"✅ Saved colorstrip with {len(order_color_map)} order colors to: {output_file}")


# ------------------- PARSE PREVIOUS iTOL ANNOTATION -------------------

def parse_previous_colors(annotation_file):
    """Parse previous iTOL TREE_COLORS to get base_family → color mapping."""
    color_dict = {}
    with open(annotation_file, 'r') as f:
        lines = f.readlines()

    current_label = None
    for line in lines:
        line = line.strip()
        if line.startswith('#'):
            current_label = line[1:].strip().split()[0]
        elif line and not line.startswith('#'):
            parts = line.split()
            if len(parts) >= 3:
                color = parts[2]
                if current_label:
                    base_label = re.sub(r'_\d+$', '', current_label)
                    if base_label not in color_dict:
                        color_dict[base_label] = color
    return color_dict

def parse_itol_annotation_color_file(annotation_file):
    mapping = {}
    current_clade = None

    with open(annotation_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("# "):
                current_clade = line[2:].strip().split('_')[0]
            elif current_clade and line:
                parts = line.split()
                if len(parts) >= 3:
                    mapping[current_clade] = parts[2]

    return mapping

# create midpoint order annotation for ring blocks..
def parse_colorstrip_file(colorstrip_path):
    color_map = {}
    with open(colorstrip_path, "r", encoding="utf-8") as f:
        in_data = False
        for line in f:
            line = line.strip()
            if line == "DATA":
                in_data = True
                continue
            if in_data and line:
                parts = line.split(",")
                if len(parts) >= 2:
                    leaf, color = parts[0], parts[1]
                    color_map[leaf] = color
    return color_map


# ------------------- ORDER COLOUR STRIP iTOL ANNOTATION -------------------

def create_itol_order_colorstrip_from_family_clade_colors(
    df,
    family_treecolors_path,
    output_file,
    darken_level=1
):
    """
    Generate an iTOL COLORSTRIP file for Orders, reusing and darkening
    the majority Family clade color from an existing TREE_COLORS file.
    """
    from collections import defaultdict

    # Step 1: Parse TREE_COLORS to get Family → color mapping
    def parse_treecolors_mrca_to_family_and_color(filepath):
        mrca_to_family = {}
        family_to_color = {}
        current_family = None
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    current_family = line[1:].strip().split()[0]
                elif line and 'clade' in line:
                    parts = line.split()
                    mrca = parts[0].replace('|', '_')
                    color = parts[2]
                    if current_family:
                        mrca_to_family[mrca] = current_family
                        if current_family not in family_to_color:
                            family_to_color[current_family] = color
        return mrca_to_family, family_to_color

    mrca_to_family, family_to_color = parse_treecolors_mrca_to_family_and_color(family_treecolors_path)

    # Step 2: Add family from MRCA to dataframe
    df = df.copy()
    df['Family_from_MRCA'] = df['Family MRCA'].map(mrca_to_family)

    # Step 3: Get majority Family per Order
    order_to_family = {}
    for order, group in df.groupby('Order'):
        families = group['Family_from_MRCA'].dropna()
        if not families.empty:
            majority_family = families.value_counts().idxmax()
            order_to_family[order] = majority_family

    # Step 4: Assign darkened colors to each Order
    order_to_color = {}
    for order, fam in order_to_family.items():
        base_color = family_to_color.get(fam, "#AAAAAA")
        dark_color = darken_color(base_color, level=darken_level)
        order_to_color[order] = dark_color

    # Step 5: Generate COLORSTRIP lines
    header = """\
    DATASET_COLORSTRIP
    SEPARATOR COMMA
    DATASET_LABEL,Order_Annotations
    COLOR,#ff0000
    COLOR_BRANCHES,0
    BORDER_WIDTH,0.5
    BORDER_COLOR,#000000
    SHOW_STRIP_LABELS,0
    DATA
    """
    lines = [textwrap.dedent(header.strip())]

    for _, row in df.iterrows():
        leaf = row['Leaves'].strip()
        order = row['Order'].strip()
        color = order_to_color.get(order)
        if color:
            lines.append(f"{leaf},{color},")
        else:
            print(f"⚠️ No color for order '{order}'")

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("\n".join(lines))

    print(f"✅ COLORSTRIP file written: {output_file}")


def parse_family_colors_and_mrca(annotation_file):
    """
    Returns:
    - family_to_color: dict of family name → color
    - mrca_to_family: dict of MRCA node → family name
    """
    family_to_color = {}
    mrca_to_family = {}
    with open(annotation_file) as f:
        current_family = None
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                current_family = line[1:].strip().split()[0]
            elif line and 'clade' in line:
                mrca = line.split()[0].replace("|", "_")
                mrca_to_family[mrca] = current_family
                if current_family not in family_to_color:
                    family_to_color[current_family] = line.split()[2]
    return family_to_color, mrca_to_family

def get_majority_family_per_order(df):
    """
    Given a DataFrame with 'Order' and 'Family' columns,
    return a dict mapping each Order → most common Family
    """
    order_to_family = {}
    for order, group in df.groupby('Order'):
        if group['Family'].isnull().all():
            majority_family = "Undefined"
        else:
            majority_family = group['Family'].value_counts().idxmax()
        order_to_family[order] = majority_family
    return order_to_family

def assign_order_colors_from_families(order_to_family, family_to_color, level=1):
    order_to_color = {}
    for order, family in order_to_family.items():
        base_color = family_to_color.get(family, '#AAAAAA')
        darker = darken_color(base_color, level=level)
        order_to_color[order] = darker
    return order_to_color


# ----------------------------------------------------------------------------

def write_order_treecolors(df, order_to_color, output_file):
    header = """\
    TREE_COLORS
    SEPARATOR SPACE
    SHOW_LEGEND 1
    LEGEND_TITLE Order Clade Colors
    DATA
    """
    lines = [textwrap.dedent(header.strip())]

    seen = set()
    for _, row in df.iterrows():
        mrca = row['Family MRCA']
        order = row['Order']
        if pd.isna(order) or mrca in seen:
            continue
        seen.add(mrca)
        color = order_to_color.get(order, "#888888")
        mrca_fmt = mrca_to_itol_format(mrca)
        lines.append(f"# {order}")
        lines.append(f"{mrca_fmt} clade {color} normal 1")

    with open(output_file, "w") as f:
        f.write("\n".join(lines))

    print(f"✅ Saved Order-level TREE_COLORS to {output_file}")


def write_itol_text_labels(midpoint_df, color_mapping, output_path, label_rank='Order'):
    header = """\
    DATASET_TEXT
    SEPARATOR COMMA
    DATASET_LABEL,Clade labels
    COLOR,#8174A0
    MARGIN,0
    SHOW_INTERNAL,0
    ALL_LABELS_ROTATION,0
    LABELS_BELOW,1
    VERTICAL_SHIFT,0
    STRAIGHT_LABELS,0
    ALIGN_TO_TREE,0
    SIZE_FACTOR,1
    EXTERNAL_LABEL_SHIFT,0
    DATA
    """

    output_lines = [textwrap.dedent(header)]
    for _, row in midpoint_df.iterrows():
        accn = row['Leaf']
        clade = row[label_rank]
        color = color_mapping.get(clade, '#000000')
        output_lines.append(f"{accn},{clade},-1,{color},bold,17,0")

    with open(output_path, 'w') as f:
        f.write("\n".join(output_lines))

    print(f"✅ iTOL TEXT label file written to: {output_path}")


def write_itol_metadata_from_reds(df, output_path, node_col='Node', red_col='RED'):
    """
    Writes an iTOL-compatible METADATA annotation file from a RED score DataFrame.

    Parameters:
    - df (pd.DataFrame): DataFrame with columns containing node names and RED values.
    - output_path (str): Path to write the METADATA file.
    - node_col (str): Column name containing node or leaf IDs (default = 'Node').
    - red_col (str): Column name containing RED values (default = 'RED').
    
    Output:
    - Writes a METADATA file usable in iTOL with RED values as numerical metadata.
    """
    with open(output_path, "w") as f:
        f.write("METADATA\n")
        f.write("SEPARATOR COMMA\n")
        f.write("FIELD_LABELS,RED_scores\n")
        f.write("DATA\n")
        for _, row in df.iterrows():
            node = row[node_col]
            red_value = row[red_col]
            f.write(f"{node},{red_value:.6f}\n")

