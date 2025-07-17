import numpy as np
import MDAnalysis as mda
import pandas as pd
import sys

## This script calculates deuterium order parameters of phopholipid tails within a certain selection of a molecule of interest using MDAnalysis

## Usage directions:  python scd_selection.py system.gro system.xtc [nth traj out of total segmented trajectories]

## The inputs include (1) the GRO file (2) the XTC file, which should have the molecule of interest centered in the box (3) the number that indicates which trajectory slice is being evaluated.
## Results are further processed with similar files used in the VMD approach; these include get-avg-sn1.sh, get-avg-sn1.sh, header-sn1.dat, header-sn2.dat

#load system
u = mda.Universe(sys.argv[1], sys.argv[2])

lipid = 'DOPC' # select lipid
n1 = 18 # sn1 tail length
n2 = 18 # sn2 tail length

#selection statement
sel = u.select_atoms("resname " + lipid, updating=True)


#---------------------------------------#
## order param for the c3/sn-1 tail ##

df = pd.DataFrame(columns = ["carbon number", "scd"])
df_n = 0
i = 2 # start from carbon 2
while i <= n1:
    cp = sel.select_atoms("name C3" + str(i), updating=True) # C3 for sn1, C2 for sn2
    if cp.n_atoms == 0:
        print("skipping " + str(i))
        continue

    hx = sel.select_atoms("name H" + str(i) + "X and not name HX", updating=True) # X/HX for sn1, R/HR for sn2
    hy = sel.select_atoms("name H" + str(i) + "Y and not name HY", updating=True) # Y/HY for sn1, S/HS for sn2
    hz = sel.select_atoms("name H" + str(i) + "Z and not name HZ", updating=True) # Z/HZ for sn1, T/HT for sn2

    cos2_theta_sum = 0.0
    nh = 0
    for ts in u.trajectory:
        nres = cp.n_atoms
        cpx = cp.positions[:,0]
        cpy = cp.positions[:,1]
        cpz = cp.positions[:,2]

        if hx.n_atoms != 0:
            hxx = hx.positions[:,0] - cpx
            hxy = hx.positions[:,1] - cpy
            hxz = hx.positions[:,2] - cpz
            for val in range(hx.n_atoms):
                norm2 = hxx[val]**2 + hxy[val]**2 + hxz[val]**2
                cos2_theta_sum = cos2_theta_sum + (hxz[val]**2/norm2)
            nh += nres

        if hy.n_atoms != 0:
            hyx = hy.positions[:,0] - cpx
            hyy = hy.positions[:,1] - cpy
            hyz = hy.positions[:,2] - cpz
            for val in range(hy.n_atoms):
                norm2 = hyx[val]**2 + hyy[val]**2 + hyz[val]**2
                cos2_theta_sum = cos2_theta_sum + (hyz[val]**2/norm2)
            nh += nres

        if hz.n_atoms != 0:
            hzx = hz.positions[:,0] - cpx
            hzy = hz.positions[:,1] - cpy
            hzz = hz.positions[:,2] - cpz
            for val in range(hz.n_atoms):
                norm2 = hxx[val]**2 + hxy[val]**2 + hxz[val]**2
                cos2_theta_sum = cos2_theta_sum + (hxz[val]**2/norm2)
            nh += nres

    scd =-1.5*cos2_theta_sum/nh + 0.5
    df.loc[df_n] = [str(i), scd]
    df_n += 1

    i += 1

np.savetxt('c3-'+sys.argv[3]+'.dat', df.values, fmt='%s %.19f', newline='\n')  #save sn1 tail information for nth trajectory slice



#------------------------------------#
## order param for the c2/sn-2 tail ##

df = pd.DataFrame(columns = ["carbon number", "scd"])
df_n = 2 # start after rows for 2R and 2S (0th and 1st rows)

# first calculate for 2R and 2S separately
cp = sel.select_atoms("name C22", updating=True) # C3 for sn1, C2 for sn2
if cp.n_atoms == 0:
    print("skipping 2")

hx = sel.select_atoms("name H2R", updating=True) # X/HX for sn1, R/HR for sn2
hy = sel.select_atoms("name H2S", updating=True) # Y/HY for sn1, S/HS for sn2

cos2_theta_sum = 0.0
cos2_theta_sumr = 0.0
cos2_theta_sums = 0.0
nh = 0
nhr = 0
nhs = 0
for ts in u.trajectory:
    nres = cp.n_atoms
    cpx = cp.positions[:,0]
    cpy = cp.positions[:,1]
    cpz = cp.positions[:,2]

    if hx.n_atoms != 0:
        hxx = hx.positions[:,0] - cpx
        hxy = hx.positions[:,1] - cpy
        hxz = hx.positions[:,2] - cpz
        for val in range(hx.n_atoms):
            norm2r = hxx[val]**2 + hxy[val]**2 + hxz[val]**2
            cos2_theta_sumr = cos2_theta_sumr + (hxz[val]**2/norm2r)
        nhr += nres

    if hy.n_atoms != 0:
        hyx = hy.positions[:,0] - cpx
        hyy = hy.positions[:,1] - cpy
        hyz = hy.positions[:,2] - cpz
        for val in range(hy.n_atoms):
            norm2s = hyx[val]**2 + hyy[val]**2 + hyz[val]**2
            cos2_theta_sums = cos2_theta_sums + (hyz[val]**2/norm2s)
        nhs += nres

scd_2r =-1.5*cos2_theta_sumr/nhr + 0.5
scd_2s =-1.5*cos2_theta_sums/nhs + 0.5
df.loc[0] = ['2R', scd_2r]
df.loc[1] = ['2S', scd_2s]		

# now calculate from carbon 3 until end
i = 3 # start from carbon 3 (carbon 2 is handled specially; see calculation for 2R and 2S above)
while i <= n2:
    cp = sel.select_atoms("name C2" + str(i), updating=True) # C3 for sn1, C2 for sn2
    if cp.n_atoms == 0:
        print("skipping " + str(i))
        continue

    hx = sel.select_atoms("name H" + str(i) + "R and not name HR", updating=True) # X/HX for sn1, R/HR for sn2
    hy = sel.select_atoms("name H" + str(i) + "S and not name HS", updating=True) # Y/HY for sn1, S/HS for sn2
    hz = sel.select_atoms("name H" + str(i) + "T and not name HT", updating=True) # Z/HZ for sn1, T/HT for sn2
    h9 = sel.select_atoms("name H" + str(i) + "1 and not name H1", updating=True) # alternate naming of H9R and H10R in some lipids

    cos2_theta_sum = 0.0
    nh = 0
    for ts in u.trajectory:
        nres = cp.n_atoms
        cpx = cp.positions[:,0]
        cpy = cp.positions[:,1]
        cpz = cp.positions[:,2]

        if hx.n_atoms != 0:
            hxx = hx.positions[:,0] - cpx
            hxy = hx.positions[:,1] - cpy
            hxz = hx.positions[:,2] - cpz
            for val in range(hx.n_atoms):
                norm2 = hxx[val]**2 + hxy[val]**2 + hxz[val]**2
                cos2_theta_sum = cos2_theta_sum + (hxz[val]**2/norm2)
            nh += nres

        if hy.n_atoms != 0:
            hyx = hy.positions[:,0] - cpx
            hyy = hy.positions[:,1] - cpy
            hyz = hy.positions[:,2] - cpz
            for val in range(hy.n_atoms):
                norm2 = hyx[val]**2 + hyy[val]**2 + hyz[val]**2
                cos2_theta_sum = cos2_theta_sum + (hyz[val]**2/norm2)
            nh += nres

        if hz.n_atoms != 0:
            hzx = hz.positions[:,0] - cpx
            hzy = hz.positions[:,1] - cpy
            hzz = hz.positions[:,2] - cpz
            for val in range(hz.n_atoms):
                norm2 = hxx[val]**2 + hxy[val]**2 + hxz[val]**2
                cos2_theta_sum = cos2_theta_sum + (hxz[val]**2/norm2)
            nh += nres

        if h9.n_atoms != 0:
            h9x = h9.positions[:,0] - cpx
            h9y = h9.positions[:,1] - cpy
            h9z = h9.positions[:,2] - cpz
            for val in range(h9.n_atoms):
                norm2 = h9x[val]**2 + h9y[val]**2 + h9z[val]**2
                cos2_theta_sum = cos2_theta_sum + (h9z[val]**2/norm2)
            nh += nres    	

    scd =-1.5*cos2_theta_sum/nh + 0.5
    df.loc[df_n] = [str(i), scd]
    df_n += 1

    i += 1

np.savetxt('c2-'+sys.argv[3]+'.dat', df.values, fmt='%s %.19f', newline='\n') #save sn2 tail information for nth trajectory slice

