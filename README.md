# analysiscodes
Some analysis codes for computer simulations

This package is designed for the post-analysis of computer simulations. Although it is initially designed for molecular dynamics (MD) simulations, in principle, it analyzes the trajectories from simulations, so it is applicable to any data consisting of atomic coordinates. This can include Monte Carlo (MC) simulations and ab-inito computations. To enable this feature, the only action is to modify the module 'dump' to read the configurations in the new format. Currently, it supports the simulation outcomes from the software [LAMMPS](https://lammps.sandia.gov) and [Hoomd-blue](http://glotzerlab.engin.umich.edu/hoomd-blue/).

This package requires some additional software for specific calculations. For example, the Voronoi tessellation is conducted by the [Voro++](http://math.lbl.gov/voro++/) package; the DCD type coordinates from Hoomd-blue are read by the [MDanalysis](https://www.mdanalysis.org) package. Therefore, these dependent packages should be installed at first. As mentioned above, these can be changed easily from the module of how to read-in the coordinates.

This package is also applicable for molecules, like polymer. But the center-of-mass or other 'fake' coordinates should be provided by projecting the atom type id to the molecule type id. The simplest way is just extracting the molecular type id for analysis, in the same way as for the atomic type.

A PDF manual has been included in the latest version!
