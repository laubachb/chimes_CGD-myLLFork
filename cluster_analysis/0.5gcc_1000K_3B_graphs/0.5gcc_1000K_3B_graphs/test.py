import numpy as np
import matplotlib.pyplot as plt
from sklearn.covariance import ledoit_wolf
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm
from sklearn.covariance import ledoit_wolf
from sklearn.decomposition import PCA

# Create an empty list to store arrays from each file
data_list = []
Covs = []

# Specify the directory
directory = "cluster_analysis/0.5gcc_1000K_3B_graphs/0.5gcc_1000K_3B_graphs/0.5gcc_1000K_3bR/"

# Loop over the range of integers from 50 to 74
for i in range(50, 75):
    # Print the current working directory
    #print("Current Working Directory:", os.getcwd())
    # Generate the filename
    filename = f'{directory}00{i}.3b_clu-r.txt'
    
    try:
        # Read the data from the file and append it to the list
        data = np.loadtxt(filename)
        #print(data)
        #print(np.shape(data[0]))
        sort_data = np.sort(data, axis=1)
        #print(sort_data[0])
        data_list.append(sort_data) # Account for graph invariance
        Wcov = ledoit_wolf(np.asarray(data_list[0]))[0]
        Covs.append(Wcov)
        print(np.shape(Covs))
        #print(np.shape(cov_image))
        #print(cov_image)
    except FileNotFoundError:
        print(f"File not found: {filename}")

# Concatenate the list of arrays along axis 0 (rows)
concat_data = np.concatenate(data_list, axis=0)

def logarithm(cov):
    d, V = np.linalg.eigh(cov)
    D = np.diag(np.log(d))
    logcov = np.dot(np.dot(V, D), V.T)
    return logcov

def sqrroot(cov):
    d, V = np.linalg.eigh(cov)
    D = np.diag(np.sqrt(d))
    sqrroot = np.dot(np.dot(V, D), V.T)
    return sqrroot

def expstep(cov,step):
    d, V = np.linalg.eigh(cov)
    D = np.diag(np.exp(d*step))
    expstep = np.dot(np.dot(V, D), V.T)
    return expstep

def mat_op(operation,d,V):
    return np.dot(V*operation(d),V.T)

#Initialize the Matrix Mean:
Covs = np.array(Covs)
print(np.shape(Covs))
geomean = np.mean(Covs, axis=0)
print(np.shape(geomean))

#Initialize the gradient descent step size and loss:
step = 1
norm_old = np.inf

#Set tolerance:
tol = 1e-8
norms = []

for n in range(10):

    #Compute the gradient
    geo_eval,geo_evec = np.linalg.eigh(geomean)
    geomean_inv_sqrt = mat_op(np.sqrt,1. / geo_eval,geo_evec)
    
    #Project matrices to tangent space and compute mean and norm:
    mats= [geomean_inv_sqrt.dot(cov).dot(geomean_inv_sqrt) for cov in Covs]
    log_mats = [logarithm(mat) for mat in mats]
    meanlog = np.mean(log_mats,axis = 0)
    norm = np.linalg.norm(meanlog)

    #Take step along identified geodesic to minimize loss:
    geomean_sqrt = sqrroot(geomean)
    geomean = geomean_sqrt.dot(expstep(meanlog,step)).dot(geomean_sqrt)

    # Update the norm and the step size
    if norm < norm_old:
        norm_old = norm

    elif norm > norm_old:
        step = step / 2.
        norm = norm_old

    if tol is not None and norm / geomean.size < tol:
        break
    
    norms.append(norm)
    
geo_eval,geo_evec = np.linalg.eigh(geomean)

geomean_inv_sqrt = mat_op(np.sqrt,1. / geo_eval,geo_evec)

def T_Project(geomean_inv_sqrt,cov):
    newmat = geomean_inv_sqrt.dot(cov).dot(geomean_inv_sqrt)
    T_cov = logarithm(newmat)
    return T_cov

T_covs = [T_Project(geomean_inv_sqrt,cov) for cov in Covs]

flat_covs = [np.ndarray.flatten(T_cov) for T_cov in T_covs]

print(np.shape(flat_covs))
pca = PCA()
X = pca.fit_transform(flat_covs)
print(X[:, 1])
# Assuming X is your transformed data using PCA
plt.scatter(X[:, 0], X[:, 1])
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('Scatter Plot of the First Two Principal Components')
plt.show()

# import tensorflow as tf
# from mpl_toolkits.mplot3d import Axes3D
# from sklearn.model_selection import train_test_split
# from keras import layers, models
# import pandas as pd
# from sklearn.preprocessing import MinMaxScaler
# from mpl_toolkits.mplot3d import Axes3D

# # Build the autoencoder
# def autoencoder(input_dim):
#     encoder = models.Sequential([
#         layers.Dense(60, activation='relu'),
#         layers.Dense(2, activation='linear') 
#     ])

#     decoder = models.Sequential([
#         layers.Dense(60, activation='relu'),
#         layers.Dense(input_dim, activation='sigmoid')
#     ])

#     model = models.Sequential([encoder, decoder])
#     model.compile(optimizer='adam', loss='mse')

#     return model

# Try Covariance stuff
# # Assuming 'data' is your numpy array of size (19000, 3)
# original_dim = concat_data.shape[1]

# # Create the autoencoder model
# autoencoder_model = autoencoder(original_dim)

# # Train the autoencoder
# autoencoder_model.fit(concat_data, concat_data, epochs=10, batch_size=32)

# # Encode the data into 3 dimensions
# encoded_data = autoencoder_model.predict(concat_data)

# # Plot the 3D scatter plot
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111)
# ax.scatter(encoded_data[:, 0], encoded_data[:, 1])
# ax.set_xlabel('Dimension 1')
# ax.set_ylabel('Dimension 2')
# ax.set_title('Autoencoder Output - 3D Scatter Plot')
# plt.show()



