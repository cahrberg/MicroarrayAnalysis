# Christian D. Ahrberg - BNTL, Sogang University - 2018
#
# Program opening two different image files. The first image file contains a reference dye for identifying wells
# The second image file contains the flourescent images. In the flourescent images the well intensities
# are measured at the positions identified from the first image. Afterwards concentrations are callculated
# using Poisson statistics.
#
###########################################################################################

# Importing required packages

import numpy as np
import cv2
import plotly.graph_objs as go
import plotly as py
import math
from operator import truediv

# ==============================================================================================================
# Defining Variables can be changed by user

Nimg        = 3     # Number of images
wellheight  = 25    # Heigth of microwells in um
wellradius  = 10    # Radius of microwells in um

# Variables for circle detection
dp1         = 1     # Ratio of accumulator, one means same resolution as input image, bigger numbers mean image is reduced 
minDist1    = 15    # Minimum distance between two circle centers in pixels
param11     = 1     # Threshold passed to Canny edge detector
param21     = 8     # Accumulatoir threshold for circles, the lower the more false circles are recognised
minRadius1  = 4     # Minimum radius of found circles in pixcels
maxRadius1  = 6     # Maximum radius of found circles in pixcels

# Defining appendixes required for opening images
name1 = "R"         # Appendix for reference dye images
name2 = "G"         # Appendix for fluorescence images
name3 = ".jpg"      # Appendix for image filetype
tag   = "-marked"   # Appendix for images in which found circles are marked

# ===============================================================================================================
# Defining constants that should not be changed by the user

Npos = 0.0 # Counter for positive wells, float so probability can be calculated
Nneg = 0.0 # Counter for negative wells, float so probability can be calculated

# ===============================================================================================================

# Function for Callculating flourescence intensity inside of cicle

def CircleIntensity(centerx,centery,radius,image, color):
    # Callculates the intensity indside of a circle
    #
    # Parameters:
    # centerx = x-coordinte of circle
    # centery = y-coordinate of circle
    # radius = radius of circle
    # image = image for callculating intensity (brightness saved as 8 bit, i.e. from 0 to 255
    # color = color to be measured (Blue [B], Green [G], or Red [R])
    #
    # Returns:
    # Intensity = average intensity value of circle in the tested image

    if color == "B":
        coval = 0
    elif color == "G":
        coval = 1
    elif color == "R":
        coval = 2
    else:
        print("Color not RGB, assuming color is green")
        coval = 1
    
    # Definging required parameters
    npixels = 0.0     # Count for pixels to find average
    brightness = 0.0  # Count for brightness of pixel

    # Creating square around circle
    #if (centerx >=2 and centery >=2 and centerx <= image.shape[0]-1
    for x in range(centerx-radius-2,centerx+radius+2):        # Varying through x
            for y in range(centery-radius-2,centery+radius+2):       # Varying through y
                 if (x <= image.shape[1]-1 and y <= image.shape[0]-1 and x>=0 and y>=0):  # Making sure coordinate in image
                     pixeldistance = math.sqrt((centerx - x)**2 + (centery - y)**2)     # Pythagoras to find radius from iterated pixcel to center of circle
                     if pixeldistance <  radius:                                        # If Pixel is in circle add to intensity callculation
                         pixel = image[y,x]
                         brightness = brightness + float(pixel[coval])/255                  # Updating total brightness
                         npixels = npixels + 1                                          # Updating total pixcel count
    if npixels == 0:
        npixels = 1 # Preventing error, division by zero
        
    Intensity = brightness / npixels                                                    # Callculating average intesnity of circle

    if Intensity == 0:      # Preventing division by zero for normalazation
        Intensity = 0.00000001 

    return Intensity

# ==============================================================================================================

# Function for plotting histograms
def Histogram(Data):

    # Making histogram according to plotly standard code
    data = [
        go.Histogram(
            x=Data, xbins=dict(start=0, size=0.005, end=1.01) 
        )
    ]

    layout = go.Layout(
        title='Histogram of flourescence intensity of wells',
        xaxis=dict(title='Flourescence intensity',
                   range=[0,1.01]),
        yaxis=dict(title='Well count'), 
        )

    fig = go.Figure(data=data, layout=layout)

    py.offline.plot(fig)

# ===============================================================================================================

# Function for esimating Concentration

def ConcCallculation(pHat,Npart,Vol):
    # Function callculating the concentration of outcome of dPCR
    #
    # Parameters:
    # pHat  =   estimated propability of positive partition
    # Npart =   total number of partitions
    # Vol   =   Volume of partitions in uL
    #
    # Returns:
    # C_est = callculated concentration in #particles/uL
    # C_low = lower confidence intervall of calculated concentration (95% z-distribution) in # particles/uL
    # C_upp = upper confidence intervall of calculated concentration (95% z-distribution) in # particles/uL
    #
    #######################################################

    # Defining constants
    zc = 1.96   # 95% confidence intervall z-distribution

    # Callculation of confidence interval on pHat
    pHat_Dev = zc * math.sqrt((pHat * (1-pHat))/Npart)  # Standard Deviation on expected result
    p_hat_low = pHat - pHat_Dev  # Lower bound of p_hat
    p_hat_upp = pHat + pHat_Dev  # Upper bound of p_hat

    # Callculating mean number of molecules per patition including 95%
    # confidence intervall
    lambda1 = -math.log(1-pHat)     # average number of molecules per division as per Poission distribution
    lambda_low = -math.log(1-p_hat_low)  # lower bound of average number of molecules per division
    lambda_upp = -math.log(1-p_hat_upp)  # upper bound of average number of molecules per division

    # Callculating concentrations in mol/uL from lambda values including
    # confidence intervalls
    C_est = lambda1 / Vol       # Esitmated concentration
    C_low = lambda_low / Vol    # Estimated lower bound of concentration
    C_upp = lambda_upp / Vol    # Estimated higher bound of concentration

    return C_est, C_low, C_upp

# ==============================================================================================================

# Main code file calling the previous defined functions

# Creating file for saving results. Deleting any file with the same filename already existing in the folder
file=open("Results.txt","w")
file.write("Results from PCR Experiments")
file.close()

# Opening file for saving data
file=open("Results.txt","a")
file.write("\n")
file.write("# Image" + "\t" + "Circle x-coord" + "\t" + "Circle y-coord" + "\t" + "Circle Radius" + "\t" + "Circle Intensity" + "\n")

#################################################################################################################
# Loop going through all the images taken, identifying circles, drawing cicles, and measuring intensities

for i in range(0,Nimg):
    
    # Print current image for debugging purposes
    print "Currently on image number %d" % (i+1)

    # Creating filenames for images to be opened
    nameRef = str(i+1) + name1 + name3              # Name of image to be opened for Ref dye
    nameReftag = str(i+1) + name1 + tag + name3     # Name of image to save marked circles for Ref dye
    nameFlu = str(i+1) + name2 + name3              # Name of image to be opened for Reporter dye
    nameFlutag = str(i+1) + name2 + tag + name3     # Name of image to save marked circles for Reporter dye


    # Converting images for detection of circles and measurement of flourescence
    # Reference dye image
    imgRef  = cv2.imread(nameRef,0)     # Opening image
    imgRef  = cv2.medianBlur(imgRef,5)  # Smothening image data
    imgRef2 = cv2.imread(nameRef)       # Image for marking circles in
    imgRef3 = cv2.imread(nameRef)       # Measuring reference dye intesity in
    # Fluorescence dye image
    imgFlu  = cv2.imread(nameFlu)       # Opening image
    imgFlu2 = cv2.imread(nameFlu)       # Image for marking circles in

    # Fitting circles using Hough transform from package cv2
    # function calls as follows: cv2.HoughCircles(image, method, dp, minDist, circles, param1, param2, minRadius, maxRadius)
    circles = cv2.HoughCircles(imgRef,cv2.HOUGH_GRADIENT,dp1,minDist1,param1=param11,param2=param21,minRadius=minRadius1,maxRadius=maxRadius1)

    # Drawing fitted circles to reference dye image
    circles = np.uint16(np.around(circles)) # Preparing data for plotting
    for j in circles[0,:]:
        # draw the outer circle
        cv2.circle(imgRef2,(j[0],j[1]),j[2],(0,255,0),2)
        # draw the center of the circle
        cv2.circle(imgRef2,(j[0],j[1]),2,(0,0,255),3)
    # Saving image with marked circles
    cv2.imwrite(nameReftag,imgRef2)

    # Drawing fitted circles to reporter dye image
    circles = np.uint16(np.around(circles)) # Preparing data for plotting
    for j in circles[0,:]:
        # draw the outer circle
        cv2.circle(imgFlu2,(j[0],j[1]),j[2],(0,255,0),2)
        # draw the center of the circle
        cv2.circle(imgFlu2,(j[0],j[1]),2,(0,0,255),3)
    # Saving image with marked circles
    cv2.imwrite(nameFlutag,imgFlu2)

    # Measuring intensity in marked cicles of reporter dye in each found circle
    intensity = []  # Empty list to save circle intensities into
    for j in circles[0,:]:
        intensity.append(CircleIntensity(j[0],j[1],j[2],imgFlu,name2))

    # Measuring intensity in marked cicles of reference dye in each found circle
    intensity_ref = []  # Empty list to save circle intensities into
    for j in circles[0,:]:
        intensity_ref.append(CircleIntensity(j[0],j[1],j[2],imgRef3,name1))

    # Normalising fluorescence intensity by dividing reporter dye well flourescence by reference dye well fluorescence
    intensity_norm = map(truediv, intensity, intensity_ref)

    
 
    # Saving data to output file
    k=0     # Running index to sync intensity list with circle list
    for j in circles[0,:]:
        file.write(str(i+1) + "\t" + str(j[0]) +'\t' + str(j[1]) + '\t' + str(j[2]) + '\t' + str(intensity[k]) +
                   '\t' + str(intensity_ref[k]) + '\t' + str(intensity_norm[k]) + "\n")
        k = k+1

file.close()

#################################################################################################################
# Taking measured normalised fluorescence intensities, plotting as Histogram, defining threshold, and callculing concentrations

# Importing data from file
ff = open("Results.txt",'r')
lines = ff.readlines()[2:]
f = []
# Dividing imported data into format for plotting, i.e. extracting only normalised flourescence
for x in lines:
    f.append(x.split('\t')[6])
ff.close()

# Mapping data, converting to float
plotdata = map(float, f)

# Plotting histogram, using Histogram method defined earlier
Histogram(plotdata)

# Defining manual threshold
Threshold = input('Please enter fluorescence threshold for positive call: ')

# Counting positive and negative partitions
for x in plotdata:
    if x > 1.0:                   # Removing unplausible counts
        Nneg = Nneg + 1
    elif (x < 1.0) & (x > Threshold):         # Wells with normalised fluorescence above threshold
        Npos = Npos + 1
    elif x <= Threshold:        # Wells with normalised fluorescence below threshold
        Nneg = Nneg + 1

# Calculating volume of well
VolWell = wellheight * math.pi * wellradius ** 2    # Volume of a well in cubic micrometer
VolWell = VolWell * math.pow(10,-9)                 # Volume of a well in uL

# Callculation of concentrations according to Poission distribution
Npart = Npos + Nneg # Total number of detected wells
pHat = Npos / Npart # Experimental determined propability of positve well
C_est, C_low, C_upp = ConcCallculation(pHat,Npart,VolWell) # Callculation of concentrations according to predifined function

# Printing results of analysis to console
print('Number of Poitive calls:                 {0:.2f}' .format(Npos))
print('Number of Negative calls:                {0:.2f}' .format(Nneg))
print('Propability of positive partition:       {0:.2f} \n' .format(pHat))
# Outputting results
print('Estimated concentration:                 {0:.3f} copies/uL' .format(C_est))
print('Lower bound of 95% Confidence intervall: {0:.3f} copies/uL' .format(C_low))
print('Upper bound of 95% Confidence intervall: {0:.3f} copies/uL' .format(C_upp))

# Saving results of analysis in file named concentrations.txt
file=open("Concentrations.txt","w")
file.write("Callculated concentrations from PCR experiment")
file.write("\n")
file.write('User defined threshold value:            {0:.2f}' .format(Threshold))
file.write("\n")
file.write('Number of Poitive calls:                 {0:.2f}' .format(Npos))
file.write("\n")
file.write('Number of Negative calls:                {0:.2f}' .format(Nneg))
file.write("\n")
file.write('Propability of positive partition:       {0:.2f} \n' .format(pHat))
file.write("\n")
file.write('Estimated concentration:                 {0:.3f} copies/uL' .format(C_est))
file.write("\n")
file.write('Lower bound of 95% Confidence intervall: {0:.3f} copies/uL' .format(C_low))
file.write("\n")
file.write('Upper bound of 95% Confidence intervall: {0:.3f} copies/uL' .format(C_upp))

file.close()



    
