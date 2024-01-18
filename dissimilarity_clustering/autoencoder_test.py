import tensorflow as tf
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.model_selection import train_test_split
from keras import layers, models
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

# Build the autoencoder
def autoencoder(input_dim):
    encoder = models.Sequential([
        layers.Dense(60, activation='relu'),
        layers.Dense(2, activation='linear')  # 3-dimensional output for visualization
    ])

    decoder = models.Sequential([
        layers.Dense(60, activation='relu'),
        layers.Dense(input_dim, activation='sigmoid')
    ])

    model = models.Sequential([encoder, decoder])
    model.compile(optimizer='adam', loss='mse')

    return model

# Specify the directory containing your files
directory_path = "all_pd"

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
        second_column = data[:60, 1].tolist()
        master_list.append(second_column)

        # Assign labels based on the condition
        label = (i // 25) % 13 + 1
        labels_list.append(label)

# Convert the master list to a NumPy array
master_array = np.array(master_list)

# Create a DataFrame from master_array
df = pd.DataFrame(data=master_array,
                  index=[i for i in range(master_array.shape[0])],
                  columns=['f'+str(i) for i in range(master_array.shape[1])])

df['label'] = labels_list  # Add labels as a new column

# Separate features (X) and target variable (y)
X = df.drop('label', axis=1)
y = df['label']

# Split the data into training and test sets
test_size = 0.2  # Adjust as needed
random_seed = 42  # Set a random seed for reproducibility
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_seed)

# Print the shapes of the training and test sets
print("X_train shape:", X_train.shape)
print("X_test shape:", X_test.shape)
print("y_train length:", len(y_train))
print("y_test length:", len(y_test))


# Create the autoencoder model
autoencoder_model = autoencoder(X_train.shape[1])

# Train the autoencoder
autoencoder_model.fit(X_train, X_train, epochs=10, batch_size=32)

# Encode the data into 3 dimensions
encoded_data = autoencoder_model.predict(X_train)

# Plot the 3D scatter plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.scatter(encoded_data[:, 0], encoded_data[:, 1], c=y_train, cmap='viridis')
ax.set_xlabel('Dimension 1')
ax.set_ylabel('Dimension 2')
ax.set_title('Autoencoder Output - 3D Scatter Plot')
plt.show()




