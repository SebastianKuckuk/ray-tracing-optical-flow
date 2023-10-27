import argparse
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Convert raw data to an image.")
parser.add_argument("--input", type=str, required=True, help="input raw data file path")
parser.add_argument("--output", type=str, required=True, help="output image file path")
parser.add_argument("--width", type=int, default=1024, help="image width")
parser.add_argument("--height", type=int, default=1024, help="image height")
parser.add_argument("--scale", type=str, default='off', help="scaling mode {off, fixed, auto, symAuto, absFixed, absAuto}")
parser.add_argument("--scaleMin", type=float, default=0, help="minimum value for scaling")
parser.add_argument("--scaleMax", type=float, default=1, help="maximum value for scaling")
args = parser.parse_args()

# Read raw binary data from the input file
with open(args.input, "rb") as file:
    raw_data = file.read()

# Calculate the number of bytes per pixel
num_color_channels = 3  # Assuming 3 channels for RGB (Red, Green, Blue)

# Convert the raw data to a NumPy array of bytes
data_array = np.frombuffer(raw_data, dtype=np.float32)

# Reshape the data array to match the image dimensions and channels
data_array = data_array.reshape((args.height, args.width, num_color_channels)).copy()

# scale if enabled
if 'off' == args.scale:
    pass  # nothing to do
elif 'fixed' == args.scale:
    # print(args.scaleMin, args.scaleMax)
    data_array = (data_array - args.scaleMin) / (args.scaleMax - args.scaleMin)
elif 'auto' == args.scale:
    min_val, max_val = data_array.min(), data_array.max()
    # print(min_val, max_val)
    scale = (max_val - min_val)
    data_array = (data_array - min_val) / (scale if scale != 0 else 1)
elif 'symAuto' == args.scale:
    extent = max(-data_array.min(), data_array.max())
    # print(extent)
    scale = 2 * extent
    data_array = (data_array + extent) / (scale if scale != 0 else 1)
elif 'absFixed' == args.scale:
    scale = args.scaleMax
    # print(scale)
    data_array = abs(data_array) / (scale if scale != 0 else 1)
elif 'absAuto' == args.scale:
    scale = abs(data_array).max()
    # print(scale)
    data_array = abs(data_array) / (scale if scale != 0 else 1)
else:
    raise Exception(f'Invalid scaling mode {args.scale}')

# Create an image from the data array
# plt.imshow(data_array)
# plt.axis('off')  # Turn off axis labels
# plt.savefig(args.output, dpi=300, bbox_inches='tight', pad_inches=0.0)  # Save the image
# plt.show()

img = Image.fromarray((data_array * 255).astype(np.uint8), mode='RGB')
img.save(args.output)
