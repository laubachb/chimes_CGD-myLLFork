<p style="text-align:center;">
    <img src="./doc/ChIMES_Github_logo-2.png" alt="" width="250"/>
</p>
<hr>

# ChIMES Cluster-Graph Fingerprinting
------------------------------------------------


*Note: This documentation is under still construction.*

The ChIMES Cluster-Graph Fingerprinting Software is a tool for generating ChIMES [1] Fingerprints which can be used for ML-IAM training set construction, active learning, adn uncertainty quantification.

* [1] [**link**](https://chemrxiv.org/engage/chemrxiv/article-details/67f0635bfa469535b9cbed7b) 1. Laubach B, Lindsey R. Cluster-Graph Fingerprinting: A Framework for Quantitative Analysis of Machine-Learned Interatomic Model Training and Simulation Data. ChemRxiv. 2025; doi:10.26434/chemrxiv-2025-vr0cs  This content is a preprint and has not been peer-reviewed.

The ChIMES Cluster-Graph Fingerprinting Software was developed at the University of Michigan, Ann Arbor in collaboration with Lawrence Livermore National Laboratory with funding from the NSF DGE 2241144.


<hr>

## Documentation
----------------

### Quick Start Guide

- Clone the Git Repository.
- Create a new directory with the location of where you would like fingerprints to be generated.
- Move the trajectory file of interest into this location (see notes on compatability).
- Copy the setup.in and run_cgd_fingerprinting.cmd files from the /examples folder into the newly created directory.
- Fill in hyperparameters within the setup.in file. Please note that everything must be declared. A description of each variable found below.
- Run the file run_cgd_fingerprint.cmd either with SBATCH or sh.

### Output Files

The software will generate four sets of files:
- *#*.xyzf        =
- *s.txt          =
- *s.txt          =
- *s.hist/*r.hist = 

### Trajectory File Compatability

Currently, the software only supports trajectory files with the following header setup:
<Number of Atoms in Frame>
NON_ORTHO <Box Dimension X> 0.0 0.0 0.0 <Box Dimension y> 0.0 0.0 0.0 <Box Dimension z>
<Atom Type> <x Position> <y Position> <z Position>

An example of this setup is described below:

64
NON_ORTHO 13.667195078 0.000000000 0.000000000 0.000000000 13.667195078 0.000000000 0.000000000 0.000000000
C 11.87156 13.16436 0.74303
C 12.03303 7.07236 12.54911
C 1.68816 12.41406 9.11085

### Hyperparameter Description

- WORKING_DIR      = Filepath for working directory (i.e. location where fingerprints will be generated)
- CGD_SRCDIR       = Source directory of fingerprinting files (/src/)
- HPC_ACCOUNT      = HPC account that will pay for work
- HPC_SYSTEM       = HPC system that software will be run on
- HPC_NODES        = Number of nodes requested for HPC system
- HPC_PPN          = Number of tasks to be invoked on each node
- HPC_WALLTIME     = Max walltime requested on HPC account
- TRAJPATH         = Filepath of trajectory file
- NFRAMES          = Number of frames to be fingerprinted in the trajectory file
- CUTOFF_2B        = 2-body cutoff (Angstroms)
- CUTOFF_3B        = 3-body cutoff (Angstroms)
- CUTOFF_4B        = 4-body cutoff (Angstroms)
- MORSE_LAMBDA     = 
- INNER_CUTOFF     = 
- NBINS_2B         = 
- NBINS_3B         = 
- NBINS_4B         = 
- JOBS_PER_BLOCK   = 
- TRANSFORMATION   = Choice of Morse or Direct transformation - dictates *.hist output.

<hr>

## Community
------------------------

Questions, discussion, and contributions (e.g. bug fixes, documentation, and extensions) are welcome. 


<hr>

## Contributing
------------------------

Contributions to the The ChIMES Cluster-Graph Fingerprinting Software should be made through a pull request, with ``develop`` as the destination branch. A test suite log file should be attached to the PR.  The `develop` branch has the latest contributions. Pull requests should target `develop`, and users who want the latest package versions, features, etc. can use `develop`.

<hr>


## Authors
----------------

The The ChIMES Cluster-Graph Fingerprinting Software was developed by Benjamin R. Laubach and Rebecca K. Lindsey.

Contributors can be found [here](https://github.com/LindseyLab-umich/chimes_CGD/graphs/contributors).

<hr>

## Citing
----------------

Please cite [1] when referencing ChIMES Cluster-Graph Fingerprints in <> a publication.

<hr>

