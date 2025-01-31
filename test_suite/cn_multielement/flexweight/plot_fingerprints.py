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
data_dict = {str(i).zfill(3): [] for i in range(1, 20)}  # e.g. '001', '002', ..., '019'

# Iterate through all files in the directory and process .hist files
hist_files = [filename for filename in os.listdir(data_dir) if filename.endswith(".hist")]

# Sort the files based on two criteria:
# 1. The leading integer (numeric) extracted from the file prefix (e.g., '001', '002', etc.)
# 2. The secondary type ('2b', '3b', '4b') extracted from the file name
hist_files.sort(key=lambda x: (int(x.split('-')[0]), x.split('.')[1]))  # First sort by prefix, then by '2b', '3b', '4b'

for filename in hist_files:
    print(filename)  # Optional: to see the files being processed
    # Extract the leading integer (e.g. '001', '002', etc.)
    prefix = filename.split('-')[0]  # '001', '002', etc.
    
    # Build the full path to the file
    file_path = os.path.join(data_dir, filename)

    # Read the second column and append it to the correct list
    if prefix in data_dict:
        data_dict[prefix].extend(read_second_column(file_path))

# Now, plot the data
plt.figure(figsize=(10, 6))

# Use the viridis colormap for coloring
cmap = plt.get_cmap("viridis")

# Get the total number of prefixes to create a color gradient
prefixes = list(data_dict.items())
num_prefixes = len(prefixes)

# Plot each prefix against one another
for idx, (prefix, fingerprint) in enumerate(prefixes):
    # Normalize the index for color mapping (from 0 to 1)
    color = cmap(idx / (num_prefixes - 1))  # Get color from colormap
    
    # Plot the data with the selected color
    plt.plot(fingerprint, label=f"Data {prefix}", color=color)

# Add labels and a legend
plt.xlabel("Index")
plt.ylabel("Value")
plt.title(f"Comparison of .hist Files (Alpha = {alpha_value})")
plt.legend(title="Leading Integer")

# Show the plot
plt.tight_layout()
plt.savefig("hist_comparison_plot.png")
