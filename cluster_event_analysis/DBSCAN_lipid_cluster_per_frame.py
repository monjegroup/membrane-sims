import sys
import MDAnalysis as mda
from sklearn.cluster import DBSCAN
import numpy as np
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt

# define the function to calculate the lipid cluster number
def analyze_cluster_counts_per_frame(universe, selection, eps_nm=1.5, min_samples=2, max_cluster_size=5, n_frames=2000):
    eps_angstrom = eps_nm * 10
    atoms = universe.select_atoms(selection)
    results = []

    for ts in universe.trajectory[-n_frames:]:

        # PBC wrap
        atoms.wrap(inplace=True)
        coords = atoms.positions

        # DBSCAN clustering
        db = DBSCAN(eps=eps_angstrom, min_samples=min_samples, metric='euclidean').fit(coords)
        labels = db.labels_

        # cluster size counting
        cluster_sizes = [np.sum(labels == i) for i in set(labels) if i != -1]
        cluster_counts = Counter(cluster_sizes)

        # obtain the record for one frame (including dimer to pentamer)
        frame_data = {'Frame': ts.frame}
        for k in range(2, max_cluster_size + 1):
            frame_data[f'Cluster_{k}'] = cluster_counts.get(k, 0)
        results.append(frame_data)

    return pd.DataFrame(results)

# set the input traj file and tpr file
tpr = r"\PATH_TO\membrane.tpr"
xtc = r"\PATH_TO\membrane.xtc"
output = r"\PATH_TO\cluster_time_series"
u = mda.Universe(tpr, xtc)

# select lipid and cutoff for lipid clusters
lipid_settings = {
    'SOPS': {"selection": "resname SOPS and name P", "eps_nm": 0.8},
    'PSM':  {"selection": "resname PSM and name P",  "eps_nm": 0.8},
    'SAPI': {"selection": "resname SAPI and name P", "eps_nm": 1.0}
}

# Perform analysis and save data
for lipid, params in lipid_settings.items():
    df = analyze_cluster_counts_per_frame(
        universe=u,
        selection=params["selection"],
        eps_nm=params["eps_nm"],
        min_samples=2,
        max_cluster_size=5,
        n_frames=2000
    )
    df.to_csv(output+f"\cluster_{lipid}_per_frame.dat", sep='\t', index=False)



