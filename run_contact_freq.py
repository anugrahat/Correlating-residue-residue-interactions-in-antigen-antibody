#!/usr/bin/env python3
"""
Contact frequency extraction for C1 WT pulling simulations (10 replicas).
Adapted from heat_count_surf.py for the C1 WT system.

For each replica, calculates per-frame residue-pair interaction frequencies
(count-based, 10 Angstrom cutoff) between antigen (resid 1-194) and
antibody (resid 195-421) surface residues (SASA > 0.5 nm^2).

Outputs per-replica files to /home/anugraha/c1_WT/analysis/repN/
"""

import numpy as np
import MDAnalysis as mda
import mdtraj as md
from tqdm import tqdm
import os
import sys

def get_surface_residues(topology_file, traj_file, start_resid, end_resid, sasa_threshold=0.5):
    """Identify surface residues using Shrake-Rupley SASA on frame 0."""
    traj = md.load_frame(traj_file, 0, top=topology_file)
    selection_str = f"resid {start_resid} to {end_resid} and protein"
    atom_indices = traj.topology.select(selection_str)

    if len(atom_indices) == 0:
        raise ValueError(f"No atoms found for residues {start_resid} to {end_resid}. Check the residue numbering.")

    traj_reduced = traj.atom_slice(atom_indices)
    sasa = md.shrake_rupley(traj_reduced, mode='residue')
    surface_residues = [residue.resSeq for residue, area in zip(traj_reduced.topology.residues, sasa[0]) if area > sasa_threshold]

    print(f"Surface residues from {start_resid} to {end_resid}: {len(surface_residues)} residues")
    print(f"  Residues: {surface_residues}")
    return surface_residues

def calculate_interaction_frequencies(traj_file, topology_file, group1_resids, group2_resids,
                                       sasa_threshold, distance_threshold, output_dir):
    """
    Calculate per-frame interaction frequencies between two groups of surface residues.

    Parameters:
        traj_file: path to .xtc trajectory
        topology_file: path to .gro topology
        group1_resids: [start, end] for antigen residues
        group2_resids: [start, end] for antibody residues
        sasa_threshold: SASA threshold in nm^2 for surface residue identification
        distance_threshold: distance cutoff in Angstroms for contacts
        output_dir: directory for output files
    """
    surface_residues_group1 = get_surface_residues(topology_file, traj_file, group1_resids[0], group1_resids[1], sasa_threshold)
    surface_residues_group2 = get_surface_residues(topology_file, traj_file, group2_resids[0], group2_resids[1], sasa_threshold)

    universe = mda.Universe(topology_file, traj_file)
    group1 = universe.select_atoms(f"resid {' '.join(map(str, surface_residues_group1))} and protein")
    group2 = universe.select_atoms(f"resid {' '.join(map(str, surface_residues_group2))} and protein")

    print(f"Number of residues in group1 (antigen): {len(group1.residues)}")
    print(f"Number of residues in group2 (antibody): {len(group2.residues)}")

    os.makedirs(output_dir, exist_ok=True)
    detail_interactions_file = os.path.join(output_dir, 'frame_residuepair_interaction_frequencies.dat')
    total_interactions_file = os.path.join(output_dir, 'frame_total_interactions.dat')

    with open(detail_interactions_file, 'w') as detail_file, \
         open(total_interactions_file, 'w') as total_file:

        detail_file.write("Frame ResiduePair InteractionFrequency\n")
        total_file.write("Frame TotalInteractions\n")

        for ts in tqdm(universe.trajectory, desc="Processing frames"):
            positions1 = group1.positions
            positions2 = group2.positions

            # Calculate pairwise distances between all atoms in group1 and group2
            distances = np.linalg.norm(positions1[:, np.newaxis] - positions2, axis=2)
            interactions = distances < distance_threshold

            total_interactions = np.sum(interactions)
            total_file.write(f"{ts.frame} {total_interactions}\n")

            # For each residue pair, count atom-atom contacts
            for i, residue1 in enumerate(group1.residues):
                for j, residue2 in enumerate(group2.residues):
                    interaction_freq = np.sum(interactions[group1.resids == residue1.resid, :][:, group2.resids == residue2.resid])
                    if interaction_freq > 0:
                        detail_file.write(f"{ts.frame} {residue1.resid}-{residue2.resid} {interaction_freq}\n")

    print(f"Data processing completed. Output saved to {output_dir}")

def main():
    # C1 WT system parameters
    # Antigen (chauA): GRO residues 1-194, Antibody (rest): GRO 195-421
    # NOTE: mdtraj 'resid' is 0-based index, not resSeq. GRO resSeq = resid + 1.
    # So antigen (resSeq 1-194) = resid 0-193, antibody (resSeq 195-421) = resid 194-420.
    group1_resids = [0, 193]    # antigen (GRO 1-194)
    group2_resids = [194, 420]  # antibody (GRO 195-421)
    sasa_threshold = 0.5        # nm^2
    distance_threshold = 10.0   # Angstroms

    base_analysis_dir = "/home/anugraha/c1_WT/analysis"

    # Define replica paths
    replica_dirs = {
        1: "/home/anugraha/c1_WT/pull",
    }
    for i in range(2, 11):
        replica_dirs[i] = f"/home/anugraha/c1_WT/pull_rep{i}"

    # Process specific replica if argument given, otherwise all
    if len(sys.argv) > 1:
        reps_to_process = [int(x) for x in sys.argv[1:]]
    else:
        reps_to_process = list(range(1, 11))

    for rep in reps_to_process:
        rep_dir = replica_dirs[rep]
        traj_file = os.path.join(rep_dir, "replica2.xtc")
        topo_file = os.path.join(rep_dir, "replica2.gro")
        output_dir = os.path.join(base_analysis_dir, f"rep{rep}")

        print(f"\n{'='*60}")
        print(f"Processing Replica {rep}")
        print(f"  Trajectory: {traj_file}")
        print(f"  Topology: {topo_file}")
        print(f"  Output: {output_dir}")
        print(f"{'='*60}")

        if not os.path.exists(traj_file):
            print(f"  WARNING: Trajectory file not found, skipping.")
            continue
        if not os.path.exists(topo_file):
            print(f"  WARNING: Topology file not found, skipping.")
            continue

        calculate_interaction_frequencies(
            traj_file, topo_file,
            group1_resids, group2_resids,
            sasa_threshold, distance_threshold,
            output_dir
        )

if __name__ == "__main__":
    main()
