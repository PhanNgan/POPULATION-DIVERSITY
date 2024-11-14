import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# Load the LOH, CNV, and scaffold lengths data
loh_data = pd.read_csv('loh_file.txt', sep='\t')
cnv_data = pd.read_csv('cnv_file.txt', sep='\t')
scaffold_lengths = pd.read_csv('scaffold_length_file.txt', sep='\t')

# Strip any leading/trailing whitespace from column names
loh_data.columns = loh_data.columns.str.strip()
cnv_data.columns = cnv_data.columns.str.strip()
scaffold_lengths.columns = scaffold_lengths.columns.str.strip()

# Set colors for each isolate
isolate_colors = {
    'Mg-Bali': 'blue', 'Mg-Borneo': 'green', 'Mg-Brazil': 'red', 'Mg-C21': 'orange',
    'Mg-C25': 'purple', 'Mg-L1': 'brown', 'Mg-L2': 'pink', 'Mg-Java2': 'cyan',
    'Mg-P': 'yellow', 'Mg-VN6': 'gray', 'Mg-VN11': 'violet', 'Mg-VN18': 'olive', 'Mg-VN27': 'magenta'
}

# Divide scaffolds into 10 groups (first 9 groups with 15 scaffolds, last one with 13 scaffolds)
scaffolds = loh_data['scaffold'].unique()
scaffold_groups = [scaffolds[i:i + 15] for i in range(0, len(scaffolds), 15)]

# Function to find overlaps and return layers for plotting
def get_layer(start, stop, layers):
    for layer, regions in enumerate(layers):
        if all(stop <= region_start or start >= region_stop for region_start, region_stop in regions):
            layers[layer].append((start, stop))
            return layer
    layers.append([(start, stop)])
    return len(layers) - 1

# Loop over each group of scaffolds to create individual plots
for group_idx, scaffold_group in enumerate(scaffold_groups, 1):
    fig, ax = plt.subplots(figsize=(10, 15))
    y_position = 0  # Start position for the first scaffold
    
    # Plot each scaffold in the group
    for scaffold in scaffold_group:
        print(f'Processing scaffold: {scaffold}')
        
        # Get LOH and CNV regions for the current scaffold
        scaffold_loh = loh_data[loh_data['scaffold'] == scaffold]
        scaffold_cnv = cnv_data[cnv_data['scaffold'] == scaffold]
        scaffold_length_row = scaffold_lengths[scaffold_lengths['scaffold'] == scaffold]
        
        # Check if scaffold length exists
        if scaffold_length_row.empty:
            print(f"No length found for scaffold: {scaffold}. Skipping.")
            continue

        # Get scaffold length
        scaffold_length = scaffold_length_row['length'].values[0]
        
        # Plot the scaffold line
        ax.plot([0, scaffold_length], [y_position, y_position], color='black', lw=2)
        
        # Label the scaffold on the left side with larger text
        ax.text(-5000, y_position, scaffold, verticalalignment='center', fontsize=12, 
                horizontalalignment='right', color='black')

        # Layers to manage overlapping regions
        layers = []

        # Plot LOH regions (multiple layers for overlaps)
        for isolate, color in isolate_colors.items():
            isolate_loh = scaffold_loh[scaffold_loh[isolate] == 1]
            for _, row in isolate_loh.iterrows():
                layer = get_layer(row['start'], row['stop'], layers)
                # Plot the LOH region as a line for this isolate
                ax.plot([row['start'], row['stop']], [y_position + layer + 1, y_position + layer + 1], color=color, lw=2)

        # Plot CNV regions (DUP and DEL, multiple layers for overlaps)
        for isolate, color in isolate_colors.items():
            for _, row in scaffold_cnv.iterrows():
                if row[isolate] == 0 or row[isolate] == 1:  # DEL (deletion)
                    # Find the layer for this region
                    layer = get_layer(row['start'], row['stop'], layers)
                    # Plot a triangle for deletion
                    triangle = plt.Polygon([[row['start'], y_position + layer + 1],
                                            [(row['start'] + row['stop']) / 2, y_position + layer + 2.5],
                                            [row['stop'], y_position + layer + 1]], color=color)
                    ax.add_patch(triangle)
                elif row[isolate] >= 3:  # DUP (duplication)
                    # Find the layer for this region
                    layer = get_layer(row['start'], row['stop'], layers)
                    # Plot an ellipse for duplication
                    ellipse = mpatches.Ellipse(((row['start'] + row['stop']) / 2, y_position + layer + 1.5),
                                               width=row['stop'] - row['start'], height=2, color=color)
                    ax.add_patch(ellipse)

        # Move to the next scaffold position
        y_position += len(layers) + 3  # Allow space for the next scaffold

    # Remove the y-axis, box, title, and legend
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_tick_params(labelsize=12, colors='black')

    # Set x-axis label
    ax.set_xlabel('Position on Scaffold', fontsize=14, color='black')

    # Save each plot to a separate file
    plt.savefig(f'loh_cnv_scaffolds_group_{group_idx}.png')

    # Show the plot (optional)
    plt.show()
