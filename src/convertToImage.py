import argparse
import numpy as np
import matplotlib.pyplot as plt

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Convert raw data to an image.")
parser.add_argument("--input", type=str, required=True, help="input raw data file path")
parser.add_argument("--output", type=str, required=True, help="output image file path")
parser.add_argument("--width", type=int, default=1024, help="image width")
parser.add_argument("--height", type=int, default=1024, help="image height")
args = parser.parse_args()

# Read raw binary data from the input file
with open(args.input, "rb") as file:
    raw_data = file.read()

# Calculate the number of bytes per pixel
num_color_channels = 3  # Assuming 3 channels for RGB (Red, Green, Blue)

# Convert the raw data to a NumPy array of bytes
data_array = np.frombuffer(raw_data, dtype=np.float32)

# Reshape the data array to match the image dimensions and channels
data_array = data_array.reshape((args.height, args.width, num_color_channels))

# Create an image from the data array
plt.imshow(data_array)
plt.axis('off')  # Turn off axis labels
plt.savefig(args.output, dpi=300, bbox_inches='tight', pad_inches=0.0)  # Save the image
# plt.show()
