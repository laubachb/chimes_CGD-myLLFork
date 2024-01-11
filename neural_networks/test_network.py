import numpy as np

# Create input tensor
# Set the dimensions of the tensor
tot_statepoints, tot_clusters, side_lengths = 100, 7, 3

# Create a random tensor
random_tensinputor = np.random.rand(tot_statepoints * tot_clusters, side_lengths)

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