# Load gromacs
module load foss gromacs/2021.5

# collect last 50ns of trajectory
# remember to change name of files to what you have/prefer
# -b and -e flags gives the time bounds in the xtc/trr file that you want to analyze in ps (here it is 450,000 - 500,000 ps)   
echo 'System' | gmx trjconv -f prod.xtc -s prod.tpr -b 400000 -e 500000 -o prod_last_100.xtc -pbc mol

# run loop for creating sections of trajectories
# remember to change cnt and cntmax to corresponding bounds
# important note: there will be 1 extra traj at the end, (e.g. traj26.xtc when I should have only 1-25). Just delete that
cnt=400000
cntmax=500000
val=1
 while [ ${cnt} -le ${cntmax} ]; do
 	let cntstep=${cnt}+10000  # 2000 represents how big each of the divided trajs will be (here it is 2000 ps (2 ns) long)
 	echo 'System' | gmx trjconv -f prod_last_100.xtc -s prod.tpr -b ${cnt} -e ${cntstep} -o traj${val}.xtc -pbc mol
 	let cnt=${cnt}+10000
 	let val=${val}+1
 done
