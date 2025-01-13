import os
import math
import pickle

def read_second_column(file_path):
    """
    Reads the second column from the given file.
    Returns the values as a list of strings.
    If the file doesn't exist, returns a list of NaNs.
    """
    values = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                # Split the line into columns
                columns = line.split()
                if len(columns) > 1:  # Ensure there is a second column
                    values.append(columns[1])
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        # Return a list of NaNs if the file doesn't exist
        return [math.nan]
    return values

def create_subdir_dictionary():
    """
    Creates a dictionary where keys are subdirectory names, and values are the concatenated
    second column values from specified files in those subdirectories.
    If a file doesn't exist, fills the corresponding value with NaNs.
    """
    # Get a list of all subdirectories in the current directory
    subdirectories = [d for d in os.listdir('.') if os.path.isdir(d)]
    
    # Dictionary to store the results
    subdir_dict = {}

    # Files to look for in each subdirectory
    target_files = [
        "training_data-training_data.2b_clu-s.hist",
        "training_data-training_data.3b_clu-s.hist",
        "training_data-training_data.4b_clu-s.hist",
    ]

    # Iterate through each subdirectory
    for subdir in subdirectories:
        subdir_path = os.path.join('.', subdir)
        collected_values = []

        # Process each target file in the subdirectory
        for file_name in target_files:
            file_path = os.path.join(subdir_path, file_name)
            collected_values.extend(read_second_column(file_path))

        # Add the collected values to the dictionary
        subdir_dict[subdir] = collected_values

    return subdir_dict

if __name__ == "__main__":
    result = create_subdir_dictionary()
    with open("silicon_fps.pkl", "wb") as pkl_file:
        pickle.dump(result, pkl_file)
    for subdir, values in result.items():
        print(f"{subdir}: {values}")
