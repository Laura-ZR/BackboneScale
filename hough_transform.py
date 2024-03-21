import numpy as np
import pydicom
from skimage.transform import hough_circle, hough_circle_peaks
from skimage.feature import canny
from skimage.filters import sobel
from skimage.morphology import closing, disk
from skimage.draw import circle_perimeter
import matplotlib.pyplot as plt

# Path to your DICOM file
dicom_file_path = "/Volumes/SILVERSD/Master Medical Imaging/1_course/team_challenge/Data/Scans volunteers/Volunteer 4/YS20220322_Scolistorm_V4_E2_MR_2/1201/DICOM/1.3.46.670589.11.20322.5.0.9332.2022032216323338013-1201-1-10ff0mp.dcm"

# Read the DICOM file
dicom_data = pydicom.dcmread(dicom_file_path)
# Extract the spacing between pixels
pixel_spacing = dicom_data.PixelSpacing

# Access pixel data
pixel_data = dicom_data.pixel_array

# Define parameters for canny edge detection
sigma = 4.0  # Standard deviation of the Gaussian filter
low_threshold = 0.1  # Lower bound threshold for hysteresis (0.1)
high_threshold = 0.3  # Upper bound threshold for hysteresis (0.3)

# Apply canny edge detection
#edges = canny(pixel_data, sigma=sigma, low_threshold=low_threshold, high_threshold=high_threshold)
edges=sobel(pixel_data)

# Define parameters for Hough Transform circle detection
radius_mm = 1.004
radius_pixels = int(radius_mm/pixel_spacing[0])
hough_radii = np.arange(radius_pixels-2, radius_pixels+3)
hough_res = hough_circle(edges, hough_radii)

# Find peaks in the Hough space
peaks = hough_circle_peaks(hough_res, hough_radii, total_num_peaks=5)  # You can adjust the number of peaks to find

# Draw circles on the original image
output_image = np.stack((pixel_data,) * 3, axis=-1)  # Convert to 3-channel image for visualization
for _, radius, row, col in zip(*peaks):
    circle_y, circle_x = circle_perimeter(row, col, radius)

    # Ensure circle coordinates fall within the image dimensions
    circle_y = np.clip(circle_y, 0, pixel_data.shape[0] - 1)
    circle_x = np.clip(circle_x, 0, pixel_data.shape[1] - 1)

    output_image[circle_y, circle_x] = (255, 0, 0)  # Red color for circles

# Display the original DICOM image with detected circles
plt.figure(figsize=(10, 8))

plt.subplot(1, 3, 1)
plt.imshow(pixel_data, cmap='gray')
plt.title('Original')
plt.axis('off')

plt.subplot(1, 3, 2)
plt.imshow(edges, cmap='gray')
plt.title('Sobel Detection')
plt.axis('off')

# plt.subplot(1, 4, 3)
# plt.imshow(closed_edges, cmap='gray')
# plt.title('Closing')
# plt.axis('off')

plt.subplot(1, 3, 3)
plt.imshow(output_image, cmap='gray')
plt.title('Detected Circles')
plt.axis('off')

plt.show()

# Now `output_image` contains the original DICOM image with detected circles drawn on it