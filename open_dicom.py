import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import numpy as np
import pydicom
import os

root = tk.Tk()
root.withdraw()

PathAux = filedialog.askdirectory(title="DICOM Images")
#Dcm_ROI = pydicom.dcmread(file_name)
#Dcm_im_ROI = np.copy(Dcm_ROI.pixel_array)
    
Dcm_im_all = []

for root,dirs,files in os.walk(PathAux):
    for i in files:
        if i.endswith(".dcm"):
            file_name=os.path.join(root,i)                  #Reads path from an image
            Dcm = pydicom.dcmread(file_name)
            Dcm_im = np.copy(Dcm.pixel_array)
            Dcm_im_all.append(Dcm_im)
        print('Name:',i)
        plt.imshow(Dcm_im, cmap = plt.cm.bone)
        plt.show()