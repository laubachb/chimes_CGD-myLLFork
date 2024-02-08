import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.decomposition import PCA
from tqdm import tqdm
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
from matplotlib.patches import Polygon
from scipy.stats import gaussian_kde


def coordinates_of_triangle_given_SSS(a, b, c):
    """a, b and c are lengths of the sides of a triangle"""

    A = (0, 0) # coordinates of vertex A

    B = (c, 0) # coordinates of vertex B

    C_x = b * float(b**2 + c**2 - a**2) / (2 * b * c)
    C_y = float(np.sqrt(b**2 - C_x**2)) # square root

    C = (float(C_x), float(C_y)) # coordinates of vertex C

    # print
    vertices = np.array([A,B,C])
    #print([A,B,C])
    return vertices


# Create an empty list to store arrays from each file
data_list = []

# Specify the directory
directory = "0.5gcc_1000k_3br/"

# Loop over the range of integers from 50 to 74
for i in range(50, 75):
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
    except FileNotFoundError:
        print(f"File not found: {filename}")

# Concatenate the list of arrays along axis 0 (rows)
concat_data = np.concatenate(data_list, axis=0)

# Print the shape of concatenated data
print("Shape:", np.shape(concat_data))

output_coords = []

for i in range(len(concat_data)):
    side_a, side_b, side_c = concat_data[i]
    output_coords.append(coordinates_of_triangle_given_SSS(side_a, side_b, side_c))

print(np.shape(np.asarray(output_coords)))
coords = np.asarray(output_coords)

# # Plot each row independently
# for i in range(len(coords)):
#     # Extract x and y coordinates for row i
#     x = coords[i, :, 0]  # Extract all x coordinates for row i
#     y = coords[i, :, 1]  # Extract all y coordinates for row i
    
#     # Append the first point to the end to close the polygon
#     x = np.append(x, x[0])
#     y = np.append(y, y[0])
    
#     # Compute kernel density estimate for the coordinates
#     xy = np.vstack([x, y])
#     z = gaussian_kde(xy)(xy)

#     # Plot the points and connect them with lines, coloring by density
#     plt.plot(x, y, marker='o', linestyle='-', color=plt.cm.viridis(z[0]))
    
# plt.show()

# Instantiate KMeans object
kmeans = KMeans(n_clusters=6)

# Fit the KMeans model to the data
kmeans.fit(concat_data)

# Get the cluster centers and labels
cluster_centers = kmeans.cluster_centers_
labels = kmeans.labels_

# # Create a 3D figure
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# # Plot the 3D scatter plot
# ax.scatter(concat_data[:, 0], concat_data[:, 1], concat_data[:, 2], c=labels)

# # Set labels and title
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
# ax.set_title('3D Scatter Plot')

# # Show the plot
# plt.show()

# # Iterate over each row in the array and apply the function
#     # Get the number of unique clusters
num_clusters = len(np.unique(labels))
start = 0
extra = 1
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(10, 10))
# # Plot clusters
for cluster_num in range(start, num_clusters):
    # Filter data points belonging to the current cluster
    x = []
    y = []
    xy = None
    z = None

    selected_data = concat_data[labels == cluster_num]
    output_coords = []
    coords = []

    for i in range(len(selected_data)):
        side_a, side_b, side_c = selected_data[i]
        output_coords.append(coordinates_of_triangle_given_SSS(side_a, side_b, side_c))

    print(np.shape(np.asarray(output_coords)))
    coords = np.asarray(output_coords)
    # Create a figure and subplots
    
    for i, ax in enumerate(axes.flatten()):
        # Plot each row independently
        for i in range(len(coords)):
            # Extract x and y coordinates for row i
            x = coords[i, :, 0]  # Extract all x coordinates for row i
            y = coords[i, :, 1]  # Extract all y coordinates for row i
            
            # Append the first point to the end to close the polygon
            x = np.append(x, x[0])
            y = np.append(y, y[0])
            
            # Compute kernel density estimate for the coordinates
            xy = np.vstack([x, y])
            z = gaussian_kde(xy)(xy)

            # Plot the points and connect them with lines, coloring by density
            ax.plot(x, y, marker='o', linestyle='-', color=plt.cm.viridis(z[0]))
            ax.set_title(f"Cluster {cluster_num}")  # Set title for each subplot
        
plt.show()
