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

class ClusterVisualizer:
    def __init__(self, directory="0.5gcc_1000K_3bS/", 
                 file_range=(50, 75), 
                 clustering_method='kmeans',
                 transformation='s',
                 cluster_dimension='3d'):
        self.directory = directory
        self.file_range = file_range
        self.data_list = []
        self.clustering_method = clustering_method
        self.transformation = transformation
        self.cluster_dimension = cluster_dimension

    def read_data(self, filename):
        try:
            data = np.loadtxt(filename)
            sort_data = np.sort(data, axis=1)
            self.data_list.append(sort_data)
        except FileNotFoundError:
            print(f"File not found: {filename}")

    def process_files(self):
        for i in range(*self.file_range):
            filename = f'{self.directory}00{i}.3b_clu-{self.transformation}.txt'
            self.read_data(filename)

    def visualize_clusters(self, clusters, concat_data):
        # Get the number of unique clusters
        num_clusters = len(np.unique(clusters))
        start = 0
        extra = 1

        if self.clustering_method == "hc":
            num_clusters += 1
            start = 1
            extra = 0

        if self.clustering_method == "dbscan":
            start = -1

        # Plot clusters
        for cluster_num in range(start, num_clusters):
            # Filter data points belonging to the current cluster
            selected_data = concat_data[clusters == cluster_num]

            # Create a 3D plot
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # Plot individual triangles with transparent grey lines
            for i in range(len(selected_data)):
                x = [selected_data[i, 0], 0, 0]
                y = [0, selected_data[i, 1], 0]
                z = [0, 0, selected_data[i, 2]]

                # Plot points
                ax.scatter(x, y, z, marker='o', c='grey')

                # Connect the points to form a triangle with transparent grey lines
                ax.plot([x[0], x[1], x[2], x[0]], [y[0], y[1], y[2], y[0]], [z[0], z[1], z[2], z[0]], c='grey', alpha=0.1)

            # Calculate the average triangle
            average_triangle = np.mean(selected_data, axis=0)
            x_avg = [average_triangle[0], 0, 0]
            y_avg = [0, average_triangle[1], 0]
            z_avg = [0, 0, average_triangle[2]]

            print("Number of Clusters: ", len(selected_data))
            print("Cluster Description: ", average_triangle)

            # Plot the average triangle with black lines
            ax.plot(x_avg + [x_avg[0]], y_avg + [y_avg[0]], z_avg + [z_avg[0]], c='black')

            if self.transformation == "s":
                # Plot red lines along the x, y, and z axes
                ax.plot([-1, 1], [0, 0], [0, 0], c='red', linestyle='--', linewidth=2)  # X-axis
                ax.plot([0, 0], [-1, 1], [0, 0], c='red', linestyle='--', linewidth=2)  # Y-axis
                ax.plot([0, 0], [0, 0], [-1, 1], c='red', linestyle='--', linewidth=2)  # Z-axis
                # Set axis limits to range from -1 to 1
                ax.set_xlim([-1, 1])
                ax.set_ylim([-1, 1])
                ax.set_zlim([-1, 1])
            elif self.transformation == "r":
                # Set axis limits to range from -1 to 1
                ax.set_xlim([0, 6])
                ax.set_ylim([0, 6])
                ax.set_zlim([0, 6])

            # Set axis labels
            ax.set_xlabel('Side A')
            ax.set_ylabel('Side B')
            ax.set_zlabel('Side C')

            ax.set_title(f"{cluster_num+extra} Cluster")

            # Show the plot
            plt.show()
    
    def coordinates_of_triangle_given_SSS(self, a, b, c):
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
    
    def visualize_triangles(self, clusters, concat_data):
        # Iterate over each row in the array and apply the function
         # Get the number of unique clusters
        num_clusters = len(np.unique(clusters))
        start = 0
        extra = 1

        if self.clustering_method == "hc":
            num_clusters += 1
            start = 1
            extra = 0

        if self.clustering_method == "dbscan":
            start = -1

        # Plot clusters
        for cluster_num in range(start, num_clusters):
            # Filter data points belonging to the current cluster
            selected_data = concat_data[clusters == cluster_num]
            # print(np.shape(selected_data)[0])
            if (self.clustering_method == 'dbscan') and (int(np.shape(selected_data)[0]) == int(0)):
                continue
            # Create the figure and axes outside the loop
            fig, ax = plt.subplots()

            print("Number of Clusters: ", len(selected_data))
            # Iterate over each row in the array and apply the function
            for row in selected_data:
                a, b, c = row
                coords = self.coordinates_of_triangle_given_SSS(a, b, c)
                edge_color = (0.5, 0.5, 0.5, 0.5)  # RGBA tuple for transparent grey
                p = Polygon(coords, edgecolor=edge_color, facecolor='none')
                ax.add_patch(p)
            
            # Calculate the average coordinates of the valid triangles
            average_triangle = np.mean(selected_data, axis=0)
            avg_a, avg_b, avg_c = average_triangle
            avg_coords = self.coordinates_of_triangle_given_SSS(avg_a, avg_b, avg_c)    
            print("Cluster Description: ", average_triangle)

            ax.set_title(f"{cluster_num+extra} Cluster")
            p = Polygon(avg_coords, edgecolor='k', facecolor='none')

            ax.add_patch(p)

            # Set the limits for the axes
            ax.set_xlim([0, 6])
            ax.set_ylim([0, 6])

            # Show the plot
            plt.show()
        return
    
    def perform_kmeans_clustering(self, concat_data, optimal_k=3):
        kmeans = KMeans(n_clusters=optimal_k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(concat_data)
        return labels
    
    def perform_hierarchical_clustering(self, concat_data):
        linkage_matrix = linkage(concat_data, method='average')

        # Choose a threshold to cut the dendrogram and form clusters
        optimal_cutoff = self.hc_silhouette_scores(concat_data, linkage_matrix)

        # Cut the dendrogram and get cluster assignments
        labels = fcluster(linkage_matrix,t=optimal_cutoff, criterion='distance')

        dendrogram(linkage_matrix, labels=labels, orientation='top')

        optimal_k = len(np.unique(labels))

        plt.title('Hierarchical Clustering Dendrogram')
        plt.xlabel('Data Points')
        plt.ylabel('Distance')
        plt.show()

        return labels, optimal_k
    
    def perform_dbscan_clustering(slef, concat_data, eps, min_samples=5):
        dbscan = DBSCAN(eps=eps, min_samples=min_samples)
        labels = dbscan.fit_predict(concat_data)
        optimal_k = len(np.unique(labels))
        return labels, optimal_k

    def plot_raw_data(self, concat_data, labels, optimal_k=3):
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(concat_data[:, 0], concat_data[:, 1], concat_data[:, 2], c=labels, marker='o')
        ax.set_xlabel('Side A')
        ax.set_ylabel('Side B')
        ax.set_zlabel('Side C')
        ax.set_title(f'{self.clustering_method} with {optimal_k} Clusters')

        if self.transformation == "s":
            # Plot red lines along the x, y, and z axes
            ax.plot([-1, 1], [0, 0], [0, 0], c='red', linestyle='--', linewidth=2)  # X-axis
            ax.plot([0, 0], [-1, 1], [0, 0], c='red', linestyle='--', linewidth=2)  # Y-axis
            ax.plot([0, 0], [0, 0], [-1, 1], c='red', linestyle='--', linewidth=2)  # Z-axis
            # Set axis limits to range from -1 to 1
            ax.set_xlim([-1, 1])
            ax.set_ylim([-1, 1])
            ax.set_zlim([-1, 1])
        elif self.transformation == "r":
            # Set axis limits to range from -1 to 1
            ax.set_xlim([0, 6])
            ax.set_ylim([0, 6])
            ax.set_zlim([0, 6])

        plt.show()
        
    def plot_histogram(self, labels):
        unique_labels, counts = np.unique(labels, return_counts=True)
        unique_labels += 1

        # Create a bar plot
        plt.bar(unique_labels, counts, color='skyblue')

        # Add labels and title
        plt.xlabel('Labels')
        plt.ylabel('Amounts')
        plt.title(f'Histogram with Labels and Amounts - {self.clustering_method}')

        # Show the plot
        plt.show()

    def kmeans_silhouette_scores(self, concat_data):
        # Perform k-means clustering for a range of k values
        silhouette_scores = []
        k_values = range(2, 30)  # You can adjust the range of k values

        for k in tqdm(k_values):
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            clusters = kmeans.fit_predict(concat_data)
            silhouette_scores.append(silhouette_score(concat_data, clusters))

        # Find the k value with the maximum silhouette score
        optimal_k = k_values[np.argmax(silhouette_scores)]

        # Plot the Silhouette Score curve
        plt.plot(k_values, silhouette_scores, marker='o')
        plt.title('Silhouette Score for Optimal K')
        plt.xlabel('Number of Clusters (K)')
        plt.ylabel('Silhouette Score')
        plt.show()

        return optimal_k
    
    def hc_silhouette_scores(self, concat_data, linkage_matrix):
        # Vary the cutoff parameter and calculate silhouette scores
        cutoff_values = np.arange(0.1, 2.0, 0.1)  # Adjust the range of cutoff values
        silhouette_scores = []

        for cutoff in cutoff_values:
            # Cut the dendrogram
            labels = fcluster(linkage_matrix, t=cutoff, criterion='distance')

            # Calculate silhouette score
            silhouette_score_value = silhouette_score(concat_data, labels)
            silhouette_scores.append(silhouette_score_value)
        # Plot the Silhouette Score curve
        plt.plot(cutoff_values, silhouette_scores, marker='o')
        plt.title('Silhouette Score for Hierarchical Clustering')
        plt.xlabel('Cutoff Parameter')
        plt.ylabel('Silhouette Score')
        plt.show()
        
        # Retrieve the corresponding cutoff value
        optimal_cutoff = cutoff_values[np.argmax(silhouette_scores)]

        # Print or use the optimal_cutoff as needed
        print(f"Optimal Cutoff Value: {optimal_cutoff}")

        return optimal_cutoff

    def dbscan_silhouette_scores(self, concat_data):
        # Set the range of eps values
        eps_values = np.arange(0.01, 0.6, 0.02)

        # Collect silhouette scores for each eps value
        silhouette_scores = []

        # Find the optimal eps value
        optimal_eps = None
        max_silhouette_score = -1

        for eps in tqdm(eps_values):
            dbscan = DBSCAN(eps=eps, min_samples=5)
            labels = dbscan.fit_predict(concat_data)

            # Check for noise points
            if len(set(labels)) > 1:
                silhouette_score_value = silhouette_score(concat_data, labels)
                silhouette_scores.append(silhouette_score_value)

                # Update optimal_eps if the silhouette score is higher
                if silhouette_score_value > max_silhouette_score:
                    max_silhouette_score = silhouette_score_value
                    optimal_eps = eps
            else:
                silhouette_scores.append(0)

        # Plot the results
        plt.plot(eps_values, silhouette_scores, marker='o')
        plt.title('Silhouette Score for Optimal Epsilon')
        plt.xlabel('Epsilon (eps)')
        plt.ylabel('Silhouette Score')
        plt.show()

        # Print the optimal eps value
        print(f"Optimal Epsilon (eps): {optimal_eps} with Silhouette Score: {max_silhouette_score}")

        return optimal_eps

    def perform_silhouette_scores(self, concat_data):
        if self.clustering_method == 'kmeans':
            hyperparameters = self.kmeans_silhouette_scores(concat_data)
        elif self.clustering_method == 'dbscan':
            hyperparameters = self.dbscan_silhouette_scores(concat_data)

        return hyperparameters
    
    def run(self):
        self.process_files()
        concat_data = np.concatenate(self.data_list, axis=0)
        print("Shape:", np.shape(concat_data))
        
        if self.clustering_method == 'kmeans':
            # Perform k-means clustering
            optimal_k = self.perform_silhouette_scores(concat_data)
            labels = self.perform_kmeans_clustering(concat_data, optimal_k)

        elif self.clustering_method == 'hc':
            # Perform hierarchical clustering
            labels, optimal_k = self.perform_hierarchical_clustering(concat_data)

        elif self.clustering_method == 'dbscan':
            epsilon = self. perform_silhouette_scores(concat_data)
            labels, optimal_k = self.perform_dbscan_clustering(concat_data, epsilon)
            
        # Visualize clusters
        if self.cluster_dimension == '3d':
            self.visualize_clusters(labels, concat_data)
        elif self.cluster_dimension == '2d':
            self.visualize_triangles(labels, concat_data)
        # Plot histogram
        self.plot_histogram(labels)

        # Plot raw data with labels
        self.plot_raw_data(concat_data, labels, optimal_k)

# Usage of the class
# cluster_visualizer = ClusterVisualizer(directory="0.5gcc_1000K_3bR/",
#                                        clustering_method="dbscan",
#                                        transformation='r',
#                                        cluster_dimension='2d')
# cluster_visualizer.run()
cluster_visualizer = ClusterVisualizer(directory="/Users/blaubach/chimes_CGD-myLLFork/cluster_analysis/1.0gcc_2000K_3B_graphs/1.0gcc_2000K_3B_graphs/1.0gcc_2000k_3bS/",
                                       file_range=(75, 100),
                                       clustering_method="hc",
                                       transformation='s',
                                       cluster_dimension='3d')
cluster_visualizer.run()