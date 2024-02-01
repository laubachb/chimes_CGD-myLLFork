import numpy as np
import pickle
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
import os
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from keras.optimizers import Adam, SGD, RMSprop, Adadelta, Adagrad, Adamax, Nadam, Ftrl
import matplotlib.pyplot as plt
from keras.metrics import sparse_categorical_accuracy

df = pd.read_excel("statepoint_conditions.xlsx")
df = df.drop(df.columns[0], axis=1)
phase_to_int = {"Diamond": 1,
                "Graphite": 2,
                "Liquid": 3}
df["Phase"] = df["Phase"].map(phase_to_int)
print(df.head())

# Repeat each row in the DataFrame 25 times
input_array = np.repeat(df.values, 25, axis=0)
print(np.shape(input_array))

with open("2b_all_pd", 'rb') as pickle_file:
    pd_2b = pickle.load(pickle_file)

with open("3b_all_pd", 'rb') as pickle_file:
    pd_3b = pickle.load(pickle_file)

with open("4b_all_pd", 'rb') as pickle_file:
    pd_4b = pickle.load(pickle_file)

# Combine the arrays along the second axis (axis=1)
all_array = np.concatenate((pd_2b, pd_3b, pd_4b), axis=2)
print(np.shape(all_array))
output_array = all_array.reshape(-1,all_array.shape[2])
print(np.shape(output_array))

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(input_array, output_array, test_size=0.2, random_state=42)

# Standardize the features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Build the neural network model
opt = Adam()
nn = Sequential()
nn.add(Dense(64, input_dim=X_train_scaled.shape[1], activation='softsign'))
nn.add(Dense(32, activation='relu'))
nn.add(Dense(180, activation='linear'))
nn.compile(loss='mse', optimizer=opt)

# Train the model
history = nn.fit(X_train_scaled, y_train, epochs=25, batch_size=32, validation_split=0.20)

# Make predictions on the test set
y_pred = nn.predict(X_test_scaled)

# Choose the index of the sample you want to check
index_to_check = 30

# Retrieve predicted and true values for the selected sample
predicted_value = y_pred[index_to_check]
true_value = y_test[index_to_check]

# # Compare the predicted and true values
# print("Predicted value:", predicted_value)
# print("True value:", true_value)

plt.plot(predicted_value, c='blue', label="Predicted")
plt.plot(true_value, c='red', label="True")
plt.legend()
plt.show()
