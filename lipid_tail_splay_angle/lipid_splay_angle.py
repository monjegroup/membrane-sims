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

# Select the P atom and N atom in the POPC lipid within top leaflet.
# if you hope to calculate other lipid, please update the following lipid name and names of end C atoms (e.g. C218 and C316):
p1 = u.select_atoms("resname POPC and name P")
c2 = u.select_atoms("resname POPC and name C2")
c214 = u.select_atoms("resname POPC and name C218")
c314 = u.select_atoms("resname POPC and name C316")


# define split angle
def split(x,y,z):
   xy = y.centroid() -x.centroid()
   xz = z.centroid() -x.centroid()
   split = np.arccos(np.dot(xy,xz)/((np.linalg.norm(xy))*(np.linalg.norm(xz))))
   return np.rad2deg(split)

# print(h1.n_atoms)
nselc = c2.n_atoms
print(nselc)

top = []
p_top = []

# Loop over the trajectory, please update the block of trajectory you hope to calculate, here I show last 2000 frames as an example.
for ts in u.trajectory[-2000:]:
   for ix in range(0,nselc):
       split_de = split(c2.split('residue')[ix],c214.split('residue')[ix],c314.split('residue')[ix])
       p_top.append(p1.split('residue')[ix].center_of_mass())
       top.append(split_de)
print(len(p_top))
print(p_top[0])
data = np.zeros((len(p_top),3))
for ix in range(len(p_top)):
    data[ix][0:2] = p_top[ix][0:2]
    data[ix][2] = top[ix]

np.savetxt('popc_splay_last100ns.dat', data, delimiter='\t', fmt='%2.4f', newline='\n')

