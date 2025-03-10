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

2. 
