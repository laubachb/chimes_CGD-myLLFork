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

    all_data = []

    # Process the chosen files
    for file_name in files_to_process:
        file_path = os.path.join(path, file_name)
        print(file_path)
        with open(file_path, 'r') as file:
            lines = file.readlines()
            print(f"Processing {mode} files - Normalized second column entries of {file_name}:")
            second_column_entries = []
            for line in lines[:60]:
                columns = line.split()
                if len(columns) >= 2:
                    second_column_entry = float(columns[1])
                    second_column_entries.append(second_column_entry)
            all_data.append(second_column_entries)
    print(np.shape(all_data))
    return all_data

# Example usage:
#os.chdir('..')
current_directory = os.getcwd()
print("Current Directory:", current_directory)

path_3b = "/all_pd_2b"  # Replace this with your actual path
path = current_directory + path_3b
print(path)

# Choose a specific chunk (e.g., Chunk_0) and mode ('equilibrium' or 'all')
# specific_chunk_index = 0
mode = "all"  # Change this to 'all' if needed
all_avgs_equilibrium = []
all_data = []

for specific_chunk_index in range(12):
    data = process_files(path, specific_chunk_index, mode)
    all_data.append(data)
    print(np.shape(data))
    avg_data = np.mean(data, axis = 0)

    # Print the result
    print("Average Array (shape 1x60):")
    avg_data = avg_data.reshape(-1,1)
    print(np.shape(avg_data))
    all_avgs_equilibrium.append(avg_data)

print("All data:" ,np.shape(all_data))
print("Avg. data:", np.shape(all_avgs_equilibrium))

# Write the object to a pickle file
pickle_filename = "2b_all_pd"
with open(pickle_filename, 'wb') as pickle_file:
    pickle.dump(all_data, pickle_file)

# Write the object to a pickle file
pickle_filename = "2b_avg_pd"
with open(pickle_filename, 'wb') as pickle_file:
    pickle.dump(all_avgs_equilibrium, pickle_file)

# Create the array
array_length = 300
increment_interval = 25

# Calculate the number of intervals
num_intervals = array_length // increment_interval

# Create the array using numpy's repeat and arange functions
labels = np.repeat(np.arange(num_intervals), increment_interval)[:array_length]

# Print the result
print(labels)

# Write the object to a pickle file
pickle_filename = "labels_pd"
with open(pickle_filename, 'wb') as pickle_file:
    pickle.dump(labels, pickle_file)

