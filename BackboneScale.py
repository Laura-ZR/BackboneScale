"""
BackboneScale 


BackboneScale is a software tool designed to accurately measure the spinal length 
from medical imaging data, specifically DICOM files. Developed with a focus on 
ease of use and precision, BackboneScale is intended to aid in the diagnosis 
and treatment of scoliosis. 

As such, the software outputs three different spinal lengths: anterior, posterior, and central. 
The anterior and posterior lengths are calculated from manual point placement. 
The central spinal length is obtained from the manual anterior and posterior points,
as well as by performing a semi-automatic spinal segmentation.

The GUI provides buttons to navigate through slices, select a DICOM folder, 
save segmentations, plot 3D points of segmented regions, 
and calculate  distances between user-placed and automatically calculated points.

Segmentations can be stored as nrrd files and a interactive 3D plot can be generated
to display the automatic points from the vertebra segmentations

Dependencies:
- Python 3.x
- Required Libraries: os, numpy, pydicom, tkinter, matplotlib, skimage, scipy, nrrd

Authors: [Jelle van der Pas, Koen Vat and Laura Zavala Rucio]
Date: [04-2024]

"""

import os
import numpy as np
import pydicom
import tkinter as tk
from tkinter import filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from skimage.segmentation import watershed
from skimage.filters import threshold_otsu
from scipy.spatial.distance import euclidean
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math
import nrrd

class BackboneScale:
    def __init__(self, master):
        """
        Initialize the DICOMApp with a tkinter master window.
        
        Parameters:
        master (tk.Tk): The tkinter master window.
        """
        self.master = master
        master.title("BackboneScale Application")
        
        # Initialize variables
        self.folder_path = None
        self.slices = []
        self.segmentations = []
        self.centers = []
        self.left_points = []
        self.right_points = []
        self.current_slice_index = 0
        self.seed_point = None
        self.total_spinal_length = 0.0
        
        #create GUI to visualize the slices and perform actions
        self.frame = tk.Frame(master)
        self.frame.pack(side=tk.BOTTOM)

        self.prev_button = tk.Button(self.frame, text="<< Prev", command=self.prev_slice)
        self.prev_button.pack(side='left')

        self.next_button = tk.Button(self.frame, text="Next >>", command=self.next_slice)
        self.next_button.pack(side='left')

        self.select_button = tk.Button(self.frame, text="Select Folder", command=self.select_folder)
        self.select_button.pack(side='right')

        self.save_button = tk.Button(self.frame, text="Save Segmentations", command=self.save_segmentations)
        self.save_button.pack(side='right')

        self.plot_button = tk.Button(self.frame, text="Plot 3D Points", command=self.plot_points)
        self.plot_button.pack(side='right')

        self.fig = Figure(figsize=(7, 7))
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        self.canvas.get_tk_widget().pack()

        self.canvas.mpl_connect('button_press_event', self.on_button_press)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)

        self.slice_label = tk.Label(master, text="")
        self.slice_label.pack(side=tk.TOP)
        
        self.calculate_manual_button = tk.Button(master, text="Calculate Manual Distances", command=self.calculate_manual_distances)
        self.calculate_manual_button.pack(side=tk.BOTTOM)

    def select_folder(self):
        """
        Open a file dialog to select a DICOM folder, then load the DICOM series.
        """
        self.folder_path = filedialog.askdirectory()
        if self.folder_path:
            self.load_dicom_series(self.folder_path)
            self.show_current_slice()

    def load_dicom_series(self, folder_path):
        """
        Load a DICOM series from a given folder path.

        Parameters:
        folder_path (str): Path to the folder containing DICOM files.
        """
        dicom_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path)]
        self.slices = [pydicom.dcmread(dicom_file) for dicom_file in dicom_files]
        self.slices.sort(key=lambda x: float(x.ImagePositionPatient[2]))
        first_slice = self.slices[0]
        if hasattr(first_slice, 'SliceThickness'):
            self.slice_thickness = float(first_slice.SliceThickness)
        elif hasattr(first_slice, 'PixelSpacing'):
            pixel_spacing = first_slice.PixelSpacing
            self.slice_thickness = float(pixel_spacing[0])

    def show_current_slice(self):
        """
        Display the current DICOM slice on the matplotlib canvas.
        """
        if self.slices:
            self.ax.clear()
            image = self.slices[self.current_slice_index].pixel_array
            self.render_image(image)
            self.slice_label.config(text=f"Slice: {self.current_slice_index + 1}/{len(self.slices)}")

    def render_image(self, image):
        """
        Render the given image on the matplotlib axes.

        Parameters:
        image (np.ndarray): Image array to display.
        """
        self.ax.imshow(image, cmap='gray')
        self.ax.axis('off')
        self.canvas.draw()

    def on_button_press(self, event):
        """
        Handle mouse button press events on the matplotlib canvas.

        Parameters:
        event (event): Mouse event.
        """
        # Left mouse button (red points)
        if event.button == 1:  # Left mouse button
            self.left_points.append((event.ydata, event.xdata))
            self.ax.plot(event.xdata, event.ydata, 'ro') 
            
        # Right mouse button (green points)  
        elif event.button == 3:  # Right mouse button
            self.right_points.append((event.ydata, event.xdata))
            self.ax.plot(event.xdata, event.ydata, 'go')  
            
        self.canvas.draw()

    def on_scroll(self, event):
        """
        Handle mouse scroll events on the matplotlib canvas.

        Parameters:
        event (event): Scroll event.
        """
        if event.button == 'up':
            self.seed_point = (int(event.ydata), int(event.xdata))
            self.ax.plot(event.xdata, event.ydata, 'bo')
            self.perform_segmentation()
        self.canvas.draw()

    def perform_segmentation(self):
        """
        Perform region growing segmentation based on the seed point.
        """
        if self.seed_point is not None:
            dicom_slice = self.slices[self.current_slice_index]
            image = dicom_slice.pixel_array
            segmentation = self.region_growing(image)
            if np.sum(segmentation) == 0:
                print("No voxels segmented. Click again to retry.")
                self.seed_point = None
                return
            self.segmentations.append(segmentation)
            self.centers.append(self.center)
            self.ax.imshow(segmentation, cmap='jet')
            self.ax.scatter(self.center[1], self.center[0], color='red')
            self.ax.set_title('Segmented Regions')
            self.ax.axis('off')
            self.canvas.draw()

    def threshold(self, image):
        """
        Returns a thresholded image array based on the otsu threshold value.

        Parameters:
        image (np.ndarray). Image array to be thresholded.
        """
        threshold_value = threshold_otsu(image)
        return image > threshold_value

    def window_threshold(self, image, factor=100):
        """
        Returns a windowed image array based on the selected intensity +- a factor.

        Parameters:
        image (np.ndarray). Image array to be windowed.
        factor (int) (optional, default = 100). Size of the window.
        """
        threshold_down = image[self.seed_point[0]][self.seed_point[1]] - factor
        threshold_up = image[self.seed_point[0]][self.seed_point[1]] + factor
        return np.logical_and(image > threshold_down, image < threshold_up)

    def region_growing(self, image):
        """
        Returns a watershed segmentation array based on the image and selected window point.

        Parameters:
        image (np.ndarray). Image array to be segmented.
        """
        factor = 100
        binary_image = self.window_threshold(image)
        markers = np.zeros_like(binary_image, dtype=bool)
        markers[self.seed_point[0], self.seed_point[1]] = True
        segmentation = watershed(binary_image, markers, mask=binary_image)
        self.center = np.mean(np.argwhere(segmentation), axis=0)
        i = 0
        while np.sqrt((self.center[0]-self.seed_point[0])**2 + (self.center[1]-self.seed_point[1])**2) > 40 and i < 5:
            binary_image = self.window_threshold(image, factor=factor//2)
            segmentation = watershed(binary_image, markers, mask=binary_image)
            self.center = np.mean(np.argwhere(segmentation), axis=0)
            i += 1
        if i == 5:
            print('Failed to segment region')
            return np.zeros_like(segmentation)
        return segmentation

  
    def calculate_manual_distances(self):
        """
        Calculate distances between manually placed points on subsequent slices.
        """
        if len(self.left_points) != len(self.right_points):
            print("Unequal number of points on left and right sides.")

        else:
            # Initialize variables for distances
            right_length = 0
            left_length = 0
            central_length = 0
            automatic_length = 0
            
            # Loop through each pair of points
            for i in range(1, len(self.left_points)):
                
                prev_slice = self.slices[i - 1]
                curr_slice = self.slices[i]
                
                #obtain distance between subsequent slices from DICOM file
                z_diff = abs(curr_slice.ImagePositionPatient[2] - prev_slice.ImagePositionPatient[2])
                
                self.spacing = z_diff
                
                #obtain xy coordinates of the manually placed points in the slice
                left_point = self.left_points[i]
                left_point_prev = self.left_points[i-1]
                right_point = self.right_points[i]
                right_point_prev = self.right_points[i-1]
                
                #calculate coordinates of point in the middle of the anterior and posterior points
                central_point = [(left_point[0] + right_point[0]) / 2, (left_point[1] + right_point[1]) / 2]
                central_point_prev = [(left_point_prev[0] + right_point_prev[0]) / 2, (left_point_prev[1] + right_point_prev[1]) / 2]
                automatic_center = self.centers[i]
                automatic_center_prev = self.centers[i-1]
                
                # Add z coordinate to all the points
                left_point_prev = np.array(list(left_point_prev + ((i-1)*z_diff,)))
                right_point_prev= np.array(list(right_point_prev+ ((i-1)*z_diff,)))
                automatic_center_prev= np.array(list(automatic_center_prev+ ((i-1)*z_diff,)))
                central_point_prev.append((i-1)*z_diff)
                central_point_prev = np.array(central_point_prev)
                    
                left_point = np.array(list(left_point+ ((i*z_diff),)))
                right_point= np.array(list(right_point + ((i*z_diff),)))
                automatic_center= np.array(list(automatic_center + ((i*z_diff),)))
                central_point.append(i*z_diff)
                central_point = np.array(central_point)
                
                #calculate euclidean distance between subsequent points and add to the total lengths
                curr_right_length = np.sum((right_point-right_point_prev)**2, axis=0)**(1/2)
                curr_left_length = np.sum((left_point-left_point_prev)**2, axis=0)**(1/2)
                curr_central_length = np.sum((central_point-central_point_prev)**2, axis=0)**(1/2)
                curr_automatic_length = np.sum((automatic_center-automatic_center_prev)**2, axis=0)**(1/2)                
                
                
                right_length += curr_right_length
                left_length += curr_left_length
                central_length += curr_central_length
                
                if automatic_length <= central_length*1.2:
                    automatic_length += curr_automatic_length
                    
                else:
                    automatic_length = z_diff
                    print("Too big")

            print("Left length in cm: ", left_length/10)
            print("Right length in cm: ", right_length/10)
            print("Central length in cm: ", central_length/10)
            
            print("Automatic length in cm: ", automatic_length/10)
            

    def next_slice(self):
        """
        Display the next slice in the DICOM series.
        """
        if self.slices:
            self.current_slice_index = (self.current_slice_index + 1) % len(self.slices)
            self.show_current_slice()

    def prev_slice(self):
        """
        Display the previous slice in the DICOM series.
        """
        if self.slices:
            self.current_slice_index = (self.current_slice_index - 1) % len(self.slices)
            self.show_current_slice()

    def save_segmentations(self):
        """
        Save the segmented volumes as NRRD files.
        """
        if self.segmentations:
            segmentation_volume = np.stack(self.segmentations, axis=-1)
            nrrd_path = os.path.join(self.folder_path, "segmentations.nrrd")
            nrrd.write(nrrd_path, segmentation_volume)
            
    def plot_points(self):
        """
        Plot 3D points of segmented regions.
        """
        if self.centers:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            centers = np.array(self.centers)
            z_coords = (np.arange(len(centers)) * self.spacing * 10)
            ax.scatter(centers[:, 1], centers[:, 0], z_coords, color='r')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z (mm)')
            ax.set_title('3D Plot of Central Points')
            
            # Adjust aspect ratio
            max_range = np.array([centers[:,1].max()-centers[:,1].min(), 
                                  centers[:,0].max()-centers[:,0].min(), 
                                  z_coords.max()-z_coords.min()]).max()
            
            mid_x = (centers[:,1].max()+centers[:,1].min())*0.5
            mid_y = (centers[:,0].max()+centers[:,0].min())*0.5
            mid_z = (z_coords.max()+z_coords.min())*0.5
            
            ax.set_xlim(mid_x - max_range/2, mid_x + max_range/2)
            ax.set_ylim(mid_y - max_range/2, mid_y + max_range/2)
            ax.set_zlim(mid_z - max_range/2, mid_z + max_range/2)
            
            plt.show()
            
def main():
    """
    Main function to start the application.
    """
    root = tk.Tk()
    app = BackboneScale(root)
    root.mainloop()
        

#Run the programm
    
if __name__ == "__main__":
    main()

