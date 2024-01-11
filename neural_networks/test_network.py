import numpy as np
import tensorflow as tf
import keras
from keras import layers, models
import matplotlib.pyplot as plt

# Create input tensor
# Set the dimensions of the tensor
tot_statepoints, tot_clusters, side_lengths = 100, 7, 3

# Create a random tensor
input = np.random.rand(tot_statepoints * tot_clusters, side_lengths)

# Reshape the tensor to (7x3)x100
input = input.reshape((tot_statepoints, tot_clusters, side_lengths))

# Print the shape of the tensor
print("Tensor shape:", np.shape(input))

# Create output tensor
# Set the dimensions of the tensor
y_dim = 100

# Create a random tensor
output = np.random.rand(tot_statepoints * y_dim)

# Reshape the tensor to (7x3)x100
output = output.reshape((tot_statepoints, y_dim))

# Print the shape of the tensor
print("Tensor shape:", np.shape(output))

# Build a simple neural network model
model = models.Sequential()
model.add(layers.Flatten(input_shape=(tot_clusters, side_lengths)))
model.add(layers.Dense(64, activation='relu'))
model.add(layers.Dense(y_dim, activation='linear'))  # Assuming linear activation for regression

# Compile the model
model.compile(optimizer='adam', loss='mean_squared_error')  # You can choose a different optimizer and loss function based on your problem

# Train the model and record the history
history = model.fit(input, output, epochs=100, batch_size=32, validation_split=0.2)  # You can adjust the number of epochs and batch size

# Plot RMSE after each epoch
train_rmse = np.sqrt(history.history['loss'])
val_rmse = np.sqrt(history.history['val_loss'])

epochs = range(1, len(train_rmse) + 1)

plt.plot(epochs, train_rmse, label='Training RMSE')
plt.plot(epochs, val_rmse, label='Validation RMSE')
plt.title('Training and Validation RMSE')
plt.xlabel('Epochs')
plt.ylabel('RMSE')
plt.legend()
plt.show()