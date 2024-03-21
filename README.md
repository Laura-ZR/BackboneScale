**Open DICOM:**

To run this code you need to download:

**- tkinter:** run in the command window "pip install tk"

**- pydicom:** run in the command window "pip install pydicom"

Then, WITHOUT CHANGING ANYTHING run the code directly. A user interface will pop-up asking for the DICOM folder to be opened.
Each slice will appear sequentially as Figure1 is closed.

**Hough Transform:**

The aim of this code is to identify circles. Specifically, recognize the central circle found in all vertebrae.
To do so, the image first undergoes edge detection and then, the Hough transform is applied.
The ultimate goal, if this works, would be to find the central point of the circle to place a landmark on each slice accordingly.

**Spinal length 3D slicer:**

To run this code you need to copy it on the Python interpreter of 3D slicer. 
The aim of this code is to measure the spinal lenght from manually placed landmarks in each slice with 3D slicer tool.
