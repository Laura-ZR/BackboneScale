
# BackboneScale

BackboneScale is a software tool designed to accurately measure the spinal length from medical imaging data, specifically DICOM files. Developed with a focus on ease of use and precision, BackboneScale is intended to aid in the diagnosis and treatment of scoliosis. As such, the software outputs three different lengths: anterior, posterior, and central.

## Table of contents
1. Key features
2. Installation
3. Usage
4. Limitations and guidance
5. Contributors

## Key features
- **DICOM support**: BackboneScale seamlessly integrates with DICOM files to easily implement its functionalities in this format. This feature makes the software ideal for medical imaging, as DICOM files are commonly outputted from scanners.  
- **Intuitive interface**: To allow an easier implementation among physicians, the tool includes a user-friendly interface. The program incorporates several buttons with their respective description that enables its use for all kinds of customers.  
- **Region-based segmentation**: To measure the spinal length, a preliminary segmentation of the spinal cord is required. For developing this software, the watershed algorithm was considered due to its intended use for segmenting overlapping structures. Moreover, it carries multiple advantages such as fastness, simplicity, low computational time and output of closed contours. 
- **Interactive visualization**: The interface incorporates multiple buttons that allow the user to move across the different slices included in the DICOM files. BackboneScale also integrates color coding to differentiate the anterior from the posterior landmarks placed by the user. 
- **Automated measurement**: Aside from the manual placing of landmarks, the software computes the spinal length automatically. BackboneScale computes four lengths, two of them, the anterior and posterior, directly calculated from the manual seeds. However, the remaining two lengths are computed from two central points identified by the system. One of the points is found from the middle point between the two manual seeds, while the other one is automatically identified as the centroid of the segmentations.
- **Export functionality**: BackBoneScale facilitates easy export of measurement results in standard formats. This implementation enables its easy integration with existing medical imaging workflows.

## Installation
1. Download the code; there are multiple options to do so:
  - Download a zip file from the menu on GitHub `Download zip`
  - Clone the repository 
```
git clone https://github.com/Laura-ZR/BackboneScale.git  

cd BackboneScale/
```
2. Install the pre-requirements from the `requirements.txt` file:
```
pip install -r requirements.txt
```

## Usage
1. Execute the code to open the software interface.
2. Click on the `Select Folder` button to identify the folder in which your DICOM files are stored. 
3. Place the two landmarks for the anterior and posterior lengths in one slice:
  - Anterior landmark (red): click with left button of the mouse 
  - Posterior landmark (green): click with right button of the mouse
4. Select the seed point for the watershed segmentation
  - Central landmark (not displayed): click with the scroll wheel of the mouse
5. Repeat this process in all slices. Travel through them by clicking `Next>>` and `<<Prev` buttons.
6. To obtain the spinal lengths, click on the `Calculate Manual Distances`.

## Credits
BackboneScale starts as a group project for Team Challenge course, as an initiative to help surgeons in the scoliosis treatment. Within this context, group 4 stands as authors of such tool:
- Koen Vat
- Jelle van der Pas
- Laura Zavala Rucio

