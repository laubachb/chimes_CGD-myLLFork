
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def print_coordinates_of_triangle_given_SSS(a, b, c):
    """a, b and c are lengths of the sides of a triangle"""
    
    A = (0, 0) # coordinates of vertex A

    B = (c, 0) # coordinates of vertex B

    C_x = b * (b**2 + c**2 - a**2) / (2 * b * c)
    C_y = np.sqrt(b**2 - C_x**2) # square root

    C = (float(C_x), float(C_y)) # coordinates of vertex C

    # print
    #vertices = zip(["A = ", "B = ", "C = "], [A, B, C])
    #print('\n'.join([v[0] + str(v[1]) for v in vertices]))
    vertices = np.array([A,B,C])
    #print([A,B,C])
    return vertices

# Create an empty list to store arrays from each file
data_list = []

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
    except FileNotFoundError:
        print(f"File not found: {filename}")

# Concatenate the list of arrays along axis 0 (rows)
concat_data = np.concatenate(data_list, axis=0)
# # Identify rows with NaN values
# nan_rows = np.isnan(concat_data).any(axis=1)

# # Drop rows with NaN values
# concat_data = concat_data[~nan_rows]


# Print the shape of concatenated data
print("Shape:", np.shape(concat_data))

# Create the figure and axes outside the loop
fig, ax = plt.subplots()

# Iterate over each row in the array and apply the function
for row in concat_data:
    a, b, c = row
    coords = print_coordinates_of_triangle_given_SSS(a, b, c)
    edge_color = (0.5, 0.5, 0.5, 0.5)  # RGBA tuple for transparent grey
    p = Polygon(coords, edgecolor=edge_color, facecolor='none')

    ax.add_patch(p)

# Calculate the average coordinates of the valid triangles
average_triangle = np.mean(concat_data, axis=0)
avg_a, avg_b, avg_c = average_triangle
avg_coords = print_coordinates_of_triangle_given_SSS(avg_a, avg_b, avg_c)
print(average_triangle)
print(avg_coords)
p = Polygon(coords, edgecolor='k', facecolor='none')

ax.add_patch(p)

# Set the limits for the axes
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])

# Show the plot after the loop
plt.show()