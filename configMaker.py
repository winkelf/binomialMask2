import os
import json

# Directory containing both .root and .fits files
base_dir = r"/media/fwinkel/TOSHIBA EXT/snolab_run1/new_cluster"

# Maximum number of pairs to include in JSON (set to None for no limit)
max_pairs = 200  

# Get all files in the directory
all_files = os.listdir(base_dir)

# Separate hits and mask files, only "commissioning"
hits_files = [
    f for f in all_files
    #if f.endswith(".root") and f.startswith("hits_corr_proc_") and "commissioning" and "2023" in f
    if f.endswith(".root") and f.startswith("hits_corr_proc_") and "commissioning" in f
]
mask_files = [
    f for f in all_files
    #if f.endswith(".fits") and f.startswith("mask_corr_proc_") and "commissioning" and "2023" in f
    if f.endswith(".fits") and f.startswith("mask_corr_proc_") and "commissioning" in f
]

file_pairs = []

# Make a dictionary for fast lookup of mask files by their variable part
mask_dict = {}
for m in mask_files:
    variable_part = m.replace("mask_corr_proc_", "").replace(".fits", "")
    mask_dict[variable_part] = os.path.join(base_dir, m)

# Match each hits file with its corresponding mask file
for h in hits_files:
    variable_part = h.replace("hits_corr_proc_", "").replace(".root", "")
    if variable_part in mask_dict:
        file_pairs.append({
            "fitsFile": mask_dict[variable_part],
            "rootFile": os.path.join(base_dir, h)
        })
    else:
        print(f"⚠ No matching mask file found for {h}")

# Apply maximum pairs limit
#if max_pairs is not None:
#    file_pairs = file_pairs[:max_pairs]

# Create the final JSON structure
output_data = {"filePairs": file_pairs}

# Save to JSON file
output_path = "file_pairs.json"
with open(output_path, "w") as f:
    json.dump(output_data, f, indent=2)

print(f" JSON file saved to {output_path} with {len(file_pairs)} pairs.")
