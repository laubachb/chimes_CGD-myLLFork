
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def print_coordinates_of_triangle_given_SSS(a, b, c):
    """a, b and c are lengths of the sides of a triangle"""

    A = (0, 0) # coordinates of vertex A

    B = (c, 0) # coordinates of vertex B

    C_x = b * float(b**2 + c**2 - a**2) / (2 * b * c)
    C_y = float(np.sqrt(b**2 - C_x**2)) # square root

    C = (float(C_x), float(C_y)) # coordinates of vertex C

    # print
    #vertices = zip(["A = ", "B = ", "C = "], [A, B, C])
    #print('\n'.join([v[0] + str(v[1]) for v in vertices]))
    vertices = np.array([A,B,C])
    #print([A,B,C])
    return vertices

# test -- print coordinates of triangle given sides of length 3, 4 and 5
print_coordinates_of_triangle_given_SSS(3, 4, 5)

# Create an empty list to store arrays from each file
data_list = []

# Specify the directory
directory = "cluster_analysis/0.5gcc_1000K_3B_graphs/0.5gcc_1000K_3B_graphs/0.5gcc_1000K_3bS/"

# Loop over the range of integers from 50 to 74
for i in range(50, 75):
    # Print the current working directory
    #print("Current Working Directory:", os.getcwd())
    # Generate the filename
    filename = f'{directory}00{i}.3b_clu-s.txt'
    
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

# Iterate over each row in the array and apply the function
for row in concat_data:
    a, b, c = row
    coords = print_coordinates_of_triangle_given_SSS(a, b, c)
    p = Polygon(coords, facecolor = 'k')

    fig,ax = plt.subplots()

    ax.add_patch(p)
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    plt.show()