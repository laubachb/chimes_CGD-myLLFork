import numpy as np
import pickle
import pandas as pd
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
import os
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from keras.optimizers import Adam, SGD, RMSprop, Adadelta, Adagrad, Adamax, Nadam, Ftrl
import matplotlib.pyplot as plt
from keras.metrics import sparse_categorical_accuracy


# os.chdir("..")

with open("2b_all_pd_equilibrium", 'rb') as pickle_file:
    all_2b_pd = pickle.load(pickle_file)
    all_2b_pd = np.vstack(all_2b_pd)

with open("3b_all_pd_equilibrium", 'rb') as pickle_file:
    all_3b_pd = pickle.load(pickle_file)
    all_3b_pd = np.vstack(all_3b_pd)

with open("4b_all_pd_equilibrium", 'rb') as pickle_file:
    all_4b_pd = pickle.load(pickle_file)
    all_4b_pd = np.vstack(all_4b_pd)

with open("labels_pd_equilibrium", 'rb') as pickle_file:
    labels = pickle.load(pickle_file)

stacked_pd = np.vstack((np.vstack((all_2b_pd,all_3b_pd)),all_4b_pd))
stacked_labels = np.hstack((labels, np.hstack((labels,labels))))
print(np.shape(stacked_pd))
print(np.shape(stacked_labels))

# stacked_data = np.vstack(all_pd)
# print(np.shape(stacked_data))

# Create a DataFrame
all_df = pd.DataFrame(stacked_pd)

# Create a DataFrame
all_df['label'] = stacked_labels

# Split data into features and labels
X = all_df.drop('label', axis=1)
y = all_df['label']

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Standardize the features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Build the neural network model
opt = Adam()
nn = Sequential()
nn.add(Dense(64, input_dim=X_train_scaled.shape[1], activation='softsign'))
# nn.add(Dense(32, activation='relu'))
nn.add(Dropout(0.5)) 
nn.add(Dense(32, activation='relu'))
nn.add(Dense(len(np.unique(stacked_labels)), activation='softmax'))
nn.compile(loss='sparse_categorical_crossentropy', optimizer=opt, metrics=['accuracy'])

# Train the model
history = nn.fit(X_train_scaled, y_train, epochs=250, batch_size=32, validation_split=0.20)

# Evaluate the model on the test set
loss, accuracy = nn.evaluate(X_test_scaled, y_test)
print(f"Test Accuracy: {accuracy * 100:.2f}%")

plt.plot(history.history['loss'], label='Training Loss')
plt.plot(history.history['val_loss'], label='Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Mean Squared Error')
plt.legend()
plt.show()
