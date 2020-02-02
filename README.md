# analysiscodes
Some analysis codes for computer simulations

This package is designed for the post-analysis of computer simulations. Although it is initially designed for molecular dynamics (MD) simulations, in principle, it analyzes the trajectories from simulations, so it is applicable to any data consisting of atomic coordinates. This can include Monte Carlo (MC) simulations and ab-inito computations. To enable this feature, the only action is to modify the module 'dump' to read the configurations in the new format. Currently, it supports the simulation outcomes from the software [LAMMPS](https://lammps.sandia.gov) and Hoomd-blue.
