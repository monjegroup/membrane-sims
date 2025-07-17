import numpy as np
import MDAnalysis
import sys

# USAGE *py [gro file] [traj file]

# This script run over last 50ns of trajectory as a example, modify the frame selection (line31) as needed

#################
# TRAJECTORY
#################

u = MDAnalysis.Universe(sys.argv[1], sys.argv[2])
print(len(u.trajectory))

# selection P atom of different lipids in the top leaflet
sope = u.select_atoms("(resname SOPE) and (name P) and (prop z >= 55)", updating=True)
sops = u.select_atoms("(resname SOPS) and (name P) and (prop z >= 55)", updating=True)
sapi = u.select_atoms("(resname SAPI) and (name P) and (prop z >= 55)", updating=True)
psm = u.select_atoms("(resname PSM) and (name P) and (prop z >= 55)", updating=True)


sopet = []
sopst = []
sapit = []
psmt = []

for ts in u.trajectory[-1000:]:
   for ix in sope:
       d = ix.position
       sopet.append((d[0], d[1], d[2]))
   for ix in sops:
       d = ix.position
       sopst.append((d[0], d[1], d[2]))
   for ix in sapi:
       d = ix.position
       sapit.append((d[0], d[1], d[2]))
   for ix in psm:
       d = ix.position
       psmt.append((d[0], d[1], d[2]))


sopet = np.array(sopet)
sopst = np.array(sopst)
sapit = np.array(sapit)
psmt = np.array(psmt)


# SAVE DATA
np.savetxt('sope.dat', sopet, delimiter='\t', fmt='%2.4f', newline='\n')
np.savetxt('sops.dat', sopst, delimiter='\t', fmt='%2.4f', newline='\n')
np.savetxt('sapi.dat', sapit, delimiter='\t', fmt='%2.4f', newline='\n')
np.savetxt('psm.dat', psmt, delimiter='\t', fmt='%2.4f', newline='\n')