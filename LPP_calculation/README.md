Run Trajectory for LPP Calculation
1. Generate a trajectory for LPP calculation
   Start from a well-equilibrated configuration. Use prod_mem.mdp with PMF cutoff = 2 and nstxout = 2500 (i.e., saving every 5 ps) to run a 100 ns production trajectory:
   $gmx mdrun -deffnm prod_mem

2. Prepare TPR file for 10 ns blocks
   To split the 100 ns trajectory into 10 ns segments for LPP calculation, generate a .tpr file with the following settings in gromacs_prod.mdp:
   nsteps = 5000000 (10 ns with 2 fs timestep)
   Ensure wall time < 72h (or depends on the request of your computational cluster)
   export gmx="/projects/academic/vmonje/software/gromacs-ls/gromacs-ls-2016.3/my_install/bin/gmx_LS"
   $gmx grompp -f gromacs_prod.mdp -c step7.gro -r step7.gro -p topol.top -o traj10ns.tpr -maxwarn -2

3. Split trajectory into 10 ns segments
   Use trjconv to extract each 10 ns segment from the 100 ns trajectory:
   $gmx trjconv -f prod_mem.trr -s traj10ns.tpr -o traj_01.trr -b 0 -e 10000 -center
   Repeat for all time intervals: 0–10 ns, 10–20 ns, ..., 90–100 ns.

4. Set up each segment for LPP calculation
   For each segment (e.g., traj_01.trr), create a separate folder (traj_01/) and copy the LPP script:
   mkdir traj_01
   cp /projects/academic/vmonje/jinhuili/scripts/gromacs_ls/gromacs_ls.csh traj_01/
   In the gromacs_ls.csh, the main command to analyze the LPP is:
   srun -n 1 $gmx mdrun -s ../traj20ns.tpr -rerun ../traj_20.trr -localsgrid 0.1 -lsgridx 1 -lsgridy 1 -lsmindihang 0.0005
   Submit the job inside each folder.

Notes:
    Use .trr files instead of .xtc, because .trr preserves full atomic information (positions, velocities, and forces), which is necessary for LPP calculation.
    Each LPP job will generate a .trr output. Keeping each job in a separate folder prevents overwriting and simulation crashes.
    You will obtain the LPP for each 10 ns interval. In total, submit 10 jobs and average the results using Tensortools.
	GROMACS_LS doesn't support parallel computation. For more details, please check the gromacs_ls link: https://www.mdstress.org/index.php/gromacs-ls/

Run LPP Analysis
1. Set up Tensortools
   export tensor="/projects/academic/vmonje/software/gromacs-ls/gromacs-4.5.5-ls-5.0/my_install/bin/tensortools"

2. Average the LPP outputs
   $tensor -f 1.dat0,2.dat0,...,10.dat0 -o avg.dat0     # Average over 10 segments
   $tensor -f avg.dat0 --prof z -o stress.dat --gf 2     # Calculate stress profile along z-axis

3. Post-process LPP profile
   Compute lateral pressure profile using:
   awk '{printf "%.1f %.3f\n", $1, -($2+$6)/2 - $10}' stress.dat > stress_fixed.dat
   Columns: $2 = Pxx, $6 = Pyy, $10 = Pzz
   Formula: LPP = –(Pxx + Pyy)/2 – Pzz

