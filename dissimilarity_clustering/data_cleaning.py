import pickle
import os.path

path_3b = "/all_pd_3b/"

current_directory = os.getcwd()
print("Current Directory:", current_directory)

path = current_directory + path_3b

# List all files in the directory
files = os.listdir(path)

# Print the list of files
print(files)


