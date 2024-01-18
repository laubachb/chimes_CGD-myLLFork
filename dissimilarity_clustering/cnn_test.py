from keras.models import Sequential
from keras.layers import Dense
from keras.optimizers import Adam
from keras.utils import to_categorical
from sklearn.model_selection import train_test_split
import numpy as np
import pandas as pd
import os
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
from sklearn import preprocessing
import tensorflow as tf
from tensorflow.keras import layers

# Specify the directory containing your files
directory_path = "../dissimilarity_clustering/all_pd_3b"

# Initialize the master list
master_list = []
labels_list = []  # List to store labels

# Iterate through each file in the directory
for i, filename in enumerate(os.listdir(directory_path)):
    if filename.endswith('.hist'):
        # Construct the full file path
        file_path = os.path.join(directory_path, filename)

        # Load the data from the file
        data = np.loadtxt(file_path)

        # Extract the second column into a list and append to the master list
        second_column = data[:60].tolist()
        master_list.append(second_column)

        # Assign labels based on the condition
        label = (i // 25) % 13 + 1
        labels_list.append(label)

# Convert the master list to a NumPy array
master_array = np.array(master_list)

print(np.shape(master_array))
print(np.shape(labels_list))



