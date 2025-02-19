import os
import numpy as np
import matplotlib.pyplot as plt

# Path to your directory containing the .hist files
data_dir = os.getcwd()
setup_file = os.path.join(data_dir, "setup.in")

# Function to extract alpha value from setup.in
def get_alpha(setup_path):
    try:
        with open(setup_path, 'r') as file:
            for line in file:
                if "alpha" in line.lower():
                    return line.split('=')[-1].strip()
    except FileNotFoundError:
        print("setup.in not found, using default alpha value")
    return "Unknown"

# Function to read second column of a .hist file
def read_second_column(file_path):
    data = np.loadtxt(file_path)
    return data[:, 1]  # Assuming second column is the one we need

alpha_value = get_alpha(setup_file)

# Dictionary to hold lists of second column values for each leading integer
data_dict = {str(i).zfill(4): [] for i in range(1, 20)}  # e.g. '001', '002', ..., '019'

# Iterate through all files in the directory and process .hist files
hist_files = [filename for filename in os.listdir(data_dir) if filename.endswith(".hist")]

# Sort the files based on two criteria:
# 1. The leading integer (numeric) extracted from the file prefix (e.g., '001', '002', etc.)
# 2. The secondary type ('2b', '3b', '4b') extracted from the file name
# First, sort the hist_files by prefix (as int) and then by the identifier (e.g., '2b', '3b', etc.)
hist_files.sort(key=lambda x: (int(x.split('-')[0]), x.split('.')[1]))
print("Hist. Files (sorted):", hist_files)

# Rotate the list so that the file with prefix "0000" comes first.
start_index = 0
for i, filename in enumerate(hist_files):
    if filename.split('-')[0] == "0000":
        start_index = i
        break
hist_files = hist_files[start_index:] + hist_files[:start_index]

print("Hist. Files (rotated to start at 0000):", hist_files)

# Process each file
for filename in hist_files:
    print("Processing file:", filename)
    # Extract the leading integer (e.g. '0000', '0001', etc.)
    prefix = filename.split('-')[0]
    
    # Build the full path to the file
    file_path = os.path.join(data_dir, filename)

    # Read the second column and append it to the correct list in data_dict.
    # (Assumes that data_dict already has keys for these prefixes.)
    if prefix in data_dict:
        data_dict[prefix].extend(read_second_column(file_path))

# Filter out any prefixes whose lists are empty.
data_dict = {k: v for k, v in data_dict.items() if v}

# Now, plot the data.
plt.figure(figsize=(10, 6))
print("Data to be plotted:", data_dict)

# Use the viridis colormap for coloring.
cmap = plt.get_cmap("viridis")

# Create a list of (prefix, fingerprint) pairs from data_dict.
prefixes = list(data_dict.items())
num_prefixes = len(prefixes)
print("Number of nonempty prefixes:", num_prefixes)

# Plot each prefix's data with a unique color.
for idx, (prefix, fingerprint) in enumerate(prefixes):
    # Normalize the index for color mapping (from 0 to 1)
    color = cmap(idx / (num_prefixes - 1)) if num_prefixes > 1 else cmap(0)
    plt.plot(fingerprint, label=f"Data {prefix}", color=color)

# Add labels, a title, and a legend.
plt.xlabel("Index")
plt.ylabel("Value")
plt.title(f"Comparison of .hist Files (Alpha = {alpha_value})")
plt.legend(title="Leading Integer")

# Finalize and save the plot.
plt.tight_layout()
plt.savefig("hist_comparison_plot.png")
