import pickle
import os.path
import numpy as np

import os

def process_files(path, specific_chunk_index, mode="all"):
    # List all files in the directory
    files = os.listdir(path)

    # Ignore files that start with "03"
    filtered_files = [file for file in files if not file.startswith("03")]

    # Sort the filtered files based on file names
    sorted_files = sorted(filtered_files)

    # Separate into chunks of 25 based on file names
    chunk_size = 25
    file_chunks = [sorted_files[i:i + chunk_size] for i in range(0, len(sorted_files), chunk_size)]

    # Choose a specific chunk
    specific_chunk = file_chunks[specific_chunk_index]

    # Determine which files to consider based on the mode
    if mode == "equilibrium":
        start_index = len(specific_chunk) // 2
        files_to_process = specific_chunk[start_index:]
    elif mode == "all":
        files_to_process = specific_chunk
    else:
        raise ValueError("Invalid mode. Use 'equilibrium' or 'all'.")

    # Process the chosen files
    for file_name in files_to_process:
        file_path = os.path.join(path, file_name)
        with open(file_path, 'r') as file:
            lines = file.readlines()
            print(f"Processing {mode} files - Normalized second column entries of {file_name}:")
            second_column_entries = []
            for line in lines[:60]:
                columns = line.split()
                if len(columns) >= 2:
                    second_column_entry = float(columns[1])
                    second_column_entries.append(second_column_entry)

            # Normalize the second column entries between 0 and 1
            min_value = min(second_column_entries)
            max_value = max(second_column_entries)

            normalized_entries = [(x - min_value) / (max_value - min_value) for x in second_column_entries]

            # Print the normalized entries
            print(normalized_entries)

# Example usage:
current_directory = os.getcwd()
print("Current Directory:", current_directory)

path_3b = "/your/path/to/3b"  # Replace this with your actual path
path = os.path.join(current_directory, path_3b)

# Choose a specific chunk (e.g., Chunk_0) and mode ('equilibrium' or 'all')
specific_chunk_index = 0
mode = "equilibrium"  # Change this to 'all' if needed

process_files(path, specific_chunk_index, mode)
