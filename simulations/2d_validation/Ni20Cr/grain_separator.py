import cv2
import numpy as np
import matplotlib.pyplot as plt
import math
# import random as rng
# import itertools
import datetime
from string import Template

img_name = 'Ni20Cr_EBSD.bmp'
img = cv2.imread(img_name)

### GETTING THE SCALE OF THE IMAGE
scale_num = 100 #um
units       = 'microns'
w_int = 1
n_per_w = 5
h_width = 2*w_int/n_per_w
H = 100
W = 150
l_PF = 0.5

gray = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
ret,thresh = cv2.threshold(gray,5,255,cv2.THRESH_BINARY)

contours, hierarchy = cv2.findContours(thresh, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)

#typically the first contour is the scale THRESH_BINARY
x,y,w,h=cv2.boundingRect(contours[0])

dpp = scale_num/(w-2)    #distance per pixel



print("Image stepsizes are: ",dpp ," units/pixel")



### CROP THE IMAGE TO REMOVE SCALE BAR

#Contour detects all white areas, the one with the largest area should be the scale bar. If this doesn't work, use manual cropping
ret,thresh = cv2.threshold(gray,254,255,cv2.THRESH_BINARY)

contours, hierarchy = cv2.findContours(thresh, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)

max_area=0
max_id = 0
for (idx,cnt) in enumerate(contours):
    area = cv2.contourArea(cnt)
    if area > max_area:
        max_area = area
        max_id = idx
x,y,w,h=cv2.boundingRect(contours[max_id])

cropped = img[:y,:]

print("Image original domain size is",np.asarray(cropped.shape[:-1])*dpp)



cropped = cropped[:int(H/dpp),:int(W/dpp)]
cv2.imwrite('cropped_microstructure.png',cropped)

### SEPARATING GRAINS USING HUE

hue = cv2.cvtColor(cropped,cv2.COLOR_BGR2GRAY)#hsv[:,:,0]

h,w = hue.shape

hist,bins = np.histogram(hue.ravel(),256,[0,256])

### IDENTIFYING HUE PEAKS

peak_ids = []
for idx in range(0,hist.shape[0]):
    try:
        if ((hist[idx] > hist[idx-1]) and (hist[idx] > hist[idx+1])):
            peak_ids+=[idx]
        elif ((hist[idx] > hist[idx+1]) and idx ==0):
            peak_ids+=[idx]
        elif ((hist[idx] > hist[idx-1]) and idx == hist.shape[0]-2):
            peak_ids+=[idx]
    except:
        continue

hue_peaks = bins[peak_ids]

### USING FLOOD FILL TO IDENTIFY GRAIN NUMBERS

grains = np.zeros((h,w),np.uint16)

grn_ct = 0
GRN_list = []

for thresh in range(256):#hue_peaks:
    thresh_hue = np.asarray(255*(hue==thresh),np.uint8)
    contours, hierarchy = cv2.findContours(thresh_hue, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

    GRN_list+=[cnt for cnt in contours]

cntsSorted = sorted(GRN_list, key=lambda x: cv2.contourArea(x),reverse=True)

A = [cv2.contourArea(cnt) for cnt in cntsSorted]
D = [2*math.sqrt(a_i/math.pi) for a_i in A]
print("Average grain diameter = ", 2*math.sqrt(sum(A)/len(A)/math.pi)*dpp, " ",units," Number of grains = ",len(A) )
print("Range of diameters",min(D),max(D))

cv2.drawContours(cropped,cntsSorted,-1,[0,0,0],1)
tmp = cv2.flip(cropped,0)

for (idx,cnt) in enumerate(cntsSorted):
    cv2.fillPoly(grains,pts=[cnt],color=[idx])
    cv2.drawContours(grains, [cnt], 0, [idx], 1)

### GENERATING X,Y AND GRAIN NUMBER OUTPUT
uniqueValues, occurCount = np.unique(grains, return_counts=True)
phases = np.zeros(grains.shape)

molten_salt_region = (uniqueValues.shape[0])*np.ones((h,int( (W+10)/dpp ) ),np.uint16)

grains = np.hstack([grains,molten_salt_region])
h,w = grains.shape
phases = np.hstack([phases,np.ones(molten_salt_region.shape)])

print("Final mesh elements numbers: ",grains.shape)
print("Number of grains ",uniqueValues.shape[0])

print("Mesh size in pixels = ",h," X ",w)
print("Mesh dimensions = ",h*dpp," X ",w*dpp)
print("Mesh dpp = ",dpp)


grain_count = str(uniqueValues.shape[0]+1)
phi1        = str(0)
phi2        = str(0)
phi3        = str(0)
z           = str(0)
phase_id    = str(0)
symm_class  = str(0)
mat_name    = "Ni_5Cr" #Doesn't matter, just added for convention
time_now    = datetime.datetime.now()

err = 1e-5

x_min       = str(0 - err)
y_min       = str(0 - err)
z_min       = str(0)
x_max       = str(w*dpp + err)
y_max       = str(h*dpp + err)
z_max       = str(0)
z_step      = str(0)
z_dim       = str(0)

### GET HEADER TEMPLATE
with open('ebsd_template.txt','r') as file:
    lines = (file.readlines() )
HEADER = Template(''.join(lines))

HEADER = HEADER.substitute(date_time=time_now,units = units,phase_id=phase_id,mat_name=mat_name,symm_class=symm_class,grn_ct = grain_count,x_min=x_min,x_max=x_max,dpp=dpp,w=w,y_min=y_min,y_max=y_max,h=h,z_min=z_min,z_max=z_max,z_step=z_step,z_dim=z_dim)

with open('EBSD_IC.txt','w') as file:
    file.write(HEADER)
###HEADER

    for i in range(0,h,1):
        for j in range(0,w,1):
            x = str((j)*dpp)
            y = str((i)*dpp)
            gr_id = str(grains[i][j])
            phase_id = str(phases[i][j])
            data = [phi1,phi2,phi3,x,y,z,gr_id,phase_id,symm_class]
            line = ' '.join(data)
            file.write(line)
            file.write("\n")
