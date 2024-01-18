import numpy as np
import pickle
import pandas as pd
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
import os

os.chdir("..")

with open("2b_all_pd_equilibrium", 'rb') as pickle_file:
    all_pd = pickle.load(pickle_file)

with open("labels_pd_equilibrium", 'rb') as pickle_file:
    labels = pickle.load(pickle_file)

print(np.shape(all_pd))

stacked_data = np.vstack(all_pd)
print(np.shape(stacked_data))

# Create a DataFrame
all_df = pd.DataFrame(stacked_data)

# Create a DataFrame
all_df['label'] = labels

# Split the data into features and labels
X = all_df.drop('label', axis=1).values
y = all_df['label'].values

# Define the neural network model
model = Sequential()
model.add(Dense(32, input_shape=(60,), activation='relu'))
model.add(Dense(16, activation='relu'))
model.add(Dense(1, activation='sigmoid'))

# Compile the model
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

# Train the model
model.fit(X, y, epochs=100, batch_size=32, validation_split=0.2)

# Evaluate the model
loss, accuracy = model.evaluate(X, y)
print(f"Accuracy: {accuracy * 100:.2f}%")



