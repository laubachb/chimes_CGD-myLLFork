# Set the input file name
set datafile separator whitespace  # Columns are separated by whitespace
filename = "000-000.2b_clu-s.hist"              # Replace with your actual file name

# Set labels for the plot
set xlabel "X-axis Label"          # Replace with your desired label for the x-axis
set ylabel "Y-axis Label"          # Replace with your desired label for the y-axis
set title "Plot of Column 1 vs Column 2"

# Customize output (optional)
set terminal pngcairo              # Output format
set output "output.png"            # Name of the output file

# Plot the data
plot filename using 1:2 with linespoints title "Data"

