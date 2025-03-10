# Description of analysis codes

1. Name of the source file: cluster_analysis.f

What it does: Solid cluster identification using both a distance criterion and the order parameter described
in Keys and Glotzer, Phys. Rev. Lett. 99, 235503 (2007), but using l=12 instead of l=10,
calculation of radial density, pair distribution function, van Hove autocorrelation function and average 
coordination number.

How to compile: ifort -O3 -132 -o cluster_analysis.x cluster_analysis.f

Files needed to run cluster_analysis.x:  It needs three files:
* in-cluster.d
* the same input file as the MC_GPU.exe code, i.e. input.d
* the movie created by the MC_GPU.exe code: movie.xyz

Output files:
* Cluster-size.dat. File with 8 columns: 1st is number of frame, 2nd is the size of the largest
cluster identified with a distance criterion, 3rd-7th are the sizes of the five largest clusters 
identified with the OP based on Steinhardt order parameters, and 8th is the total number of particles
whose environment has been identified as solid with the OP regardless of whether they are part 
of a cluster or  not.
* movie-centred.xyz. Movie in which particles are assigned different labels depending on their assembled
state, and centred at the center of mass of the largest cluster in the first frame.R
* vanHove_autorcor-tot.dat. Van Hove autocorrelation function
* PDF.dat. Radial distribution function (1st column is distance and 2nd column PDF)
* histo-vec1.dat. 1st column is coordination number, 2nd column is fraction of particles
  in the solid cluster with that coordination.
* rhoR_av.dat. Radial density averaged over all trajectory (1st column is radial distance, 2nd column is radial density).
* conf-core.xyz. Trayectory containing only those particles belong to the inner core of the solid cluster
(the radius of the sphere defining the core is provided by the user in the file in-cluster.d)
* conf-core-patch.xyz. The same as conf-core.xyz, but including the positions of the patches for each particle besides 
the position of the center of mass.

2.  Name of the source file: BOOD_bond.f 

What it does: It calculates the Bond Orientational Order Diagram, by using an energy criterium to define bonds
(two particles are considered bonded if the energy is lower than -0.2epsilon).

How to compile: ifort -O3 -132 -o BOOD_bond.x BOOD_bond.f 

Files needed to run BOOD_bond.x: in-BOOD.d and conf-core-patch.xyz 

Output files:BOOD-Lambert-av.dat    It contains the Lambert projection of the BOOD. 1st column is x coordinate,
2nd column is y coordinate, 3rd column is probability density of finding a bond with that orientation.

3.  Name of the source file: lattice_structure_factor.f

What it does: It calculates the diffration pattern projectd on the plane z=0.

How to compile: ifort -O3 -132 -o lattice_structure_factor.x lattice_structure_factor.f

Files needed to run lattice_structure_factor.x: in-pattern.d and conf-core.xyz 

Output files: Sq-av.dat    It contains the diffraction pattern projecte on the z=0 plane. 1st column is x coordinate,
2nd column is y coordinate, 3rd column is the structure factor.

4.  Name of the source file: lifting3.f

What it does: It performs the lifting to 6D.

How to compile: ifort -O3 -132 -o lifting3.x lifting3.f

Files needed to run lifting3.x: in-lift.d, bonds.xyz and coords.xyz. The file bonds.xyz contains the matching between bonds 
in 3D and bonds in the 6D lattice. coords.xyz it contains the 3D structure that will be lifted. It must be oriented so 
that 2-fold rotational axes are aligned with x-, y-, and z-. Three angles can be provided in in-lift.d to perform rotations
about x-, y- and z- axes, to get the proper orientation of the structure for lifting.

Output files:  histo_dpar_dperp_pairs.dat. The 1st column is distance in parallel space and the second column is average distance 
between pairs of particles in perp space, i.e. data needed to plot phason strain. 
coords-lifted.xyz It contains the configuration but where particles are labelled according to lifting results:
A means that the particle has been correctly lifted, M are those particles with are misaligned with respect to the bonds
in bonds.xyz tile, I means that particle does not have neighbours (isolated particle), C menas that there are conflicting
assignment for that particle.
