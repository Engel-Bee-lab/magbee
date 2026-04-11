import pandas as pd
import os

all_data = []

for file in snakemake.input.stats:
    sample_data = {}
    
    sample_name = os.path.basename(os.path.dirname(file))
    sample_data["Sample"] = sample_name
    
    with open(file, "r") as f:
        for line in f:
            line = line.strip()
            
            if not line or line.startswith("All statistics"):
                continue
            
            parts = line.split()
            
            if parts[0] == "Assembly":
                sample_data["Assembly"] = parts[1]
                continue
            
            key = " ".join(parts[:-1])
            value = parts[-1]
            
            sample_data[key] = value
    
    all_data.append(sample_data)

df = pd.DataFrame(all_data)

df = df.apply(lambda col: pd.to_numeric(col, errors="ignore"))

cols = ["Sample"] + [c for c in df.columns if c != "Sample"]
df = df[cols]

df.to_csv(snakemake.output.merged, index=False)