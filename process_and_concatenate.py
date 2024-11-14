import pandas as pd
import glob

# Define the column names since the input files do not have headers
column_names = ['Chrom', 'start', 'end', 'accept_pct', 'accepts', 'supports', 'neutrals', 'length']

# Function to concatenate fragments in a BED file
def concatenate_fragments(input_file, output_file):
    # Read the data into a DataFrame without headers
    df = pd.read_csv(input_file, sep='\t', header=None, names=column_names)

    # Prepare a list to store concatenated fragments
    concatenated = []

    # Initialize the first row for processing
    current_row = df.iloc[0]

    for i in range(1, len(df)):
        next_row = df.iloc[i]

        # Check if the chromosome is the same and the gap condition holds
        if (current_row['Chrom'] == next_row['Chrom'] and 
            (next_row['start'] - current_row['end'] < 5000) and
            (current_row['accept_pct'] > 50)):

            # Update the end position of the current row
            current_row['end'] = next_row['end']

            # Average the accept_pct
            current_row['accept_pct'] = (
                (current_row['accept_pct'] + next_row['accept_pct']) / 2
            )

            # Sum other relevant fields
            current_row['accepts'] += next_row['accepts']
            current_row['supports'] += next_row['supports']
            current_row['neutrals'] += next_row['neutrals']
            current_row['length'] += next_row['length']

        else:
            concatenated.append(current_row)
            current_row = next_row

    # Append the last row after exiting the loop
    concatenated.append(current_row)

    # Create a new DataFrame for the concatenated fragments
    result_df = pd.DataFrame(concatenated)

    # Save the result to a new file
    result_df.to_csv(output_file, index=False, sep='\t', header=False)  # No header in output
    print(f"Processed {input_file} -> {output_file}")

# Loop through each *_50_loh.merged.bed file
for input_file in glob.glob('*_50_loh.merged.bed'):
    output_file = f"{input_file[:-4]}_concatenated.bed"  # Set output filename
    concatenate_fragments(input_file, output_file)
