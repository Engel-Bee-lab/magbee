import pandas as pd
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Merge QUAST stats files")
    parser.add_argument(
        "-i", "--input", nargs="+", help="List of QUAST stats files"
    )
    parser.add_argument(
        "-o", "--output", help="Output merged CSV file"
    )
    return parser.parse_args()


# Detect if running inside Snakemake
if "snakemake" in globals():
    input_files = snakemake.input.stats
    output_file = snakemake.output.merged
else:
    args = parse_args()
    input_files = args.input
    output_file = args.output

all_data = []

for file in input_files:
    sample_data = {}
    
    # Get sample name from parent directory
    sample_name = os.path.basename(os.path.dirname(file))
    sample_data["Sample"] = sample_name
    
    with open(file, "r") as f:
        for line in f:
            line = line.strip()
            
            if not line or line.startswith("All statistics"):
                continue
            
            # More robust split (handles spaces in keys)
            if " " in line:
                key, value = line.rsplit(maxsplit=1)
            else:
                continue
            
            if key == "Assembly":
                sample_data["Assembly"] = value
            else:
                sample_data[key] = value
    
    all_data.append(sample_data)

# Create DataFrame
df = pd.DataFrame(all_data)

# Convert only numeric columns
for col in df.columns:
    if col not in ["Sample", "Assembly"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

# Reorder columns (Sample first)
cols = ["Sample"] + [c for c in df.columns if c != "Sample"]
df = df[cols]

# Write output
df.to_csv(output_file, index=False)