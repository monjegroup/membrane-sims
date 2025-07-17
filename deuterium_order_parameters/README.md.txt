1. Update the lipid name and "n1" and "n2" in the scd_calculation.py;

2. Use the section_traj_v2.csh to separate the a trajectory into ten small parts. For example, if we just analyze the last 50ns trajctory, we should first get the last 50ns, and then separate it into ten 5ns trajectories.

3. conda activate py36 #activate the python environment.

4. ./run.scd to get the scd. In the run.scd, you will see you are running scd_calculation.py for each trajctory (totally 10).

5. update the length in the head-sn1.dat and head-sn2.dat. For example, POPC, sn-1: 16. sn-2 is 18.

6. ./get-avg-sn1.dat and ./get-avg-sn2.dat to get the average and std of scd for both sn-1 and sn-2 chain.
   when you do this, make sure in the line 10 of both get-ave-sn1.dat and get-avg-sn2.dat:
   "for (( i=1; i<=10; i++ ))", i will be less equal = 10 since you have ten trajectories.