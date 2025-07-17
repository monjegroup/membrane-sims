import numpy as np
import MDAnalysis
import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.pyplot import MultipleLocator

# input gro and xtc file, please correct the file name.
u = MDAnalysis.Universe('step7.gro','step7.xtc')
print(len(u.trajectory))

# Select the P atom and N atom in the SOPE lipid within top leaflet.
# if you hope to calculate other lipid, please update the following lipid name:
h1 = u.select_atoms("prop z >= 55 and name P and resname SOPE")
t1 = u.select_atoms("prop z >= 55 and name N and resname SOPE")


# Define tilt angle function
# Tilt angle in degrees
def theta(a, b):
   ba = b.centroid() -a.centroid()
   bc = [0,0,1] #unit vector in z dir
   theta = np.arccos(np.dot(ba,bc)/((np.linalg.norm(ba))*(np.linalg.norm(bc))))
   return np.rad2deg(theta)


nselc = h1.n_atoms
print("Num of P: ", nselc)
nselc_n = t1.n_atoms
print("Num of N: ", nselc_n)

top = []
p_top = []

# Loop over the trajectory, please update the block of trajectory you hope to calculate, here I show last 1000 frames as an example.
for ts in u.trajectory[-1000:]:
   for ix in range(0,nselc):
       da = theta(h1.split('residue')[ix],t1.split('residue')[ix])
       p_top.append(h1.split('residue')[ix].center_of_mass())
       top.append(da)
print(len(p_top))
print(p_top[0])
data = np.zeros((len(p_top),4))
for ix in range(len(p_top)):
    data[ix][0:3] = p_top[ix][0:3]
    data[ix][3] = top[ix]

# save data
np.savetxt('sope_headgroup_top.dat', data, delimiter='\t', fmt='%2.4f', newline='\n')

