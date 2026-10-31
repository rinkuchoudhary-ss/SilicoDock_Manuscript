# Reproducibility package for the revision of the SilicoDock manuscript

This folder holds the derived data and the analysis scripts behind Tables S6 to S11 and the end-to-end and subgroup results added in the revision.
The generated poses, benchmark inputs, CASF-2016 scoring outputs and ERK2 case study data are in the main repository deposit.
Source code and trained weights of RIGIDiff and SMiDeN are not released (proprietary software, see Code and Data Availability).

## Data (data/)
| File | Content | Used for |
|---|---|---|
| rigidiff_pose_rmsd_casf2016.csv, rigidiff_pose_rmsd_posebusters.csv | Deviation (A) of each of the 20 RIGIDiff samples from the crystallographic pose (CASF-2016: PyMOL align; PoseBusters: RDKit GetBestRMS) | Sampling power, all tables |
| smiden_per_pose_ranking_*.csv | SMiDeN score and rank of each sample | End-to-end selection (Table S7) |
| smiden_ranking_report_*.csv | Per-complex SMiDeN ranking summary (the CASF-2016 file contains three identifiers damaged by spreadsheet formatting: 1.00E+66 = 1e66, 3.00E+92 = 3e92, 3.00E+93 = 3e93) | Table S7 |
| smiden_ranking_benchmark_complex_ids.txt | The 135 complexes of the ranking benchmark used to select SMiDeN (114 are PoseBusters complexes) | Table S7 |
| rigidiff_data_splits.json | PDB identifiers of the RIGIDiff training, validation and test subsets | Table S6, S11 |
| end_to_end_selection_*.csv | Samples with SMiDeN rank and deviation | Table S7 |
| similarity_to_training_set.csv | Sequence identity, ligand Tanimoto similarity, shared scaffolds and success per complex | Table S6 |
| pocket_similarity_to_training_set.csv | Binding-site (pocket residue) identity to the nearest training pocket and success | Table S11 |
| ligand_properties_and_success.csv | Rotatable bonds, heavy atoms, rings, ring systems, macrocycle flag and success per complex | Table S8 |
| acceptance_quantities_crystal_structures.csv | Closest ligand to protein distance, closest internal distance and largest bond-length difference per CASF-2016 crystal | Table S9 |
| threshold_sensitivity.json | Sampling power as a function of receptor distance and bond-length flag | Table S10 |

## Scripts (scripts/)
Run in the order of the numeric prefix. The paths at the top of each script point to the authors' local folders and must be edited to the location of the benchmark inputs, the generated poses and the PDBBind 2020 structures.
Requirements: Python 3, RDKit, Biopython, NumPy, pandas, SciPy.

## Run settings
200 reverse diffusion timesteps with a cosine schedule, 20 samples per complex, pocket radius 15 A, grid guidance weight 0.3, sampling temperature 0.5, random seed 42, CPU execution.
