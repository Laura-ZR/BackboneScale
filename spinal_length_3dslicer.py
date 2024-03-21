# DISCLAIMER: To get this code running, you need to paste it in the 3D Slicer Python interpreter. 
#             Otherwise it won't work, you cannot run it in Python directly.              

## OPTION 1: FOR CALCULATING DISTANCES FROM ALL POINTS TO ALL POINTS
# Get the markup point list node
pointListNode = slicer.util.getNode("F") # change F to the name of your point list 
# Get the number of control points
numControlPoints = pointListNode.GetNumberOfControlPoints()
# Initialize a list to store distances
distances = []
# Loop through each pair of control points and measure the distance
for i in range(numControlPoints):
    for j in range(i+1, numControlPoints):
        point1 = [0, 0, 0]
        point2 = [0, 0, 0]
        # Get the position of the first control point
        pointListNode.GetNthControlPointPosition(i, point1)
        # Get the position of the second control point
        pointListNode.GetNthControlPointPosition(j, point2)
        # Calculate the Euclidean distance between the two points
        distance = ((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2 + (point1[2] - point2[2])**2) ** 0.5
        # Add the distance to the list
        distances.append(distance)
        # Print the distance
        print("Distance between Control Point {} and Control Point {}: {:.2f}".format(i+1, j+1, distance))

## OPTION 2: FOR CALCULATING DISTANCES CONSECUTIVELY
# Get the markup point list node
pointListNode = slicer.util.getNode("F") # change F to the name of your point list 
# Get the number of control points
numControlPoints = pointListNode.GetNumberOfControlPoints()
# Initialize a list to store distances
distances = []
total_distance = 0
# Loop through each pair of control points and measure the distance
for i in range(numControlPoints-1):
    point1 = [-1, -1, -1]
    point2 = [0, 0, 0]
    # Get the position of the first control point
    pointListNode.GetNthControlPointPosition(i+1, point1)
    # Get the position of the second control point
    pointListNode.GetNthControlPointPosition(i, point2)
    # Calculate the Euclidean distance between the two points
    distance = ((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2 + (point1[2] - point2[2])**2) ** 0.5
    # Calculate the total distance between all points
    total_distance += distance
    # Add the distance to the list
    distances.append(distance)
    # Print the distance
    print("Distance between Control Point {} and Control Point {}: {:.2f} (mm)".format(i+1, i+2, distance))
print("Total distance: {:.2f} (mm)".format(total_distance))


