# import cv2
import numpy as np
import matplotlib.pyplot as plt
import math
# import random as rng
# import itertools
import datetime
from string import Template
import re
# import parse

file_name = 'EBSD_IC_3.txt'#'small_test.txt'
header_size = 24

salt_dim = "Z"
salt_thickness = 10 #um
bias = 1e-6

with open(file_name) as raw_file:
    header_list = [next(raw_file) for i in range(header_size) ]
    points = raw_file.read()
    # header = " ".join(header_list[1:])
    # print (header)
H_dict = {}
for line in header_list:
    F = re.match(r"# (.*): (.*)",line)
    if F:
        key,value=F.groups()
        H_dict[key]= value
H_dict[salt_dim+"_MAX"] = str(float(H_dict[salt_dim+"_MAX"]) + salt_thickness)
N_dim_old = int(H_dict[salt_dim+"_DIM"])
N_dim_new = math.floor(( float(H_dict[salt_dim+"_MAX"]) - float(H_dict[salt_dim+"_MIN"]) )/float(H_dict[salt_dim+"_STEP"]))
H_dict[salt_dim+"_DIM"] = str( N_dim_new )

# print(H_dict)

if (salt_dim == "Z"):
    X_vals = float(H_dict["X_MIN"]) + float(H_dict["X_STEP"])*np.arange(0,int(H_dict["X_DIM"]) )
    Y_vals = float(H_dict["Y_MIN"]) + float(H_dict["Y_STEP"])*np.arange(0,int(H_dict["Y_DIM"]) ) - bias
    Z_vals = float(H_dict["Z_MIN"]) + float(H_dict["Z_STEP"])*np.arange(N_dim_old,int(H_dict["Z_DIM"]) ) - bias
    H_dict["X_MIN"] = str(float(X_vals[0] - bias) )
    H_dict["Y_MIN"] = str(float(Y_vals[0] - bias) )
    H_dict["Z_MIN"] = str(float(H_dict["Z_MIN"]) -bias )
    H_dict["X_MAX"] = str(float(X_vals[-1] + bias) )
    H_dict["Y_MAX"] = str(float(Y_vals[-1] + bias) )
    H_dict["Z_MAX"] = str(float(Z_vals[-1] + bias) )

# print(min(X_vals),max(X_vals),min(Y_vals),max(Y_vals),min(Z_vals),max(Z_vals))
G = np.vstack(np.meshgrid(X_vals,Y_vals,Z_vals)).reshape(3,-1).T
# print(G.shape, 50*50*20)

S=""

###Add new header
for key in H_dict.keys():
    S = S +("# "+str(key)+ ": "+str(H_dict[key]) + "\n#\n" )

# print(S)
####Append old points
S = S + points

####Append new points
P_id = str(0)
symmetry = H_dict["Symmetry_1"]
feature_id = str(0)
phi1 = str(0)
phi2 = str(0)
phi3 = str(0)

# print(G.shape)
for i in range(G.shape[0]):
    x,y,z = G[i]
    line = " ".join([phi1,phi2,phi3,str(x),str(y),str(z),feature_id,P_id,symmetry] )
    S+=line+"\n"

with open("salt_appended.txt",'w+') as write_file:
    write_file.write(S)



# print(S)

# template = open("ebsd_template.txt").read()[2:]
