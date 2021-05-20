# MC_GPU
GPU parallel Monte Carlo code for patchy particles with the model potential described in D. Tracey, E.G. Noya, and J.P.K. Doye, The Journal of chemical physics 151, 224506 (2019).


Instructions for compilation:

ifort -axsse4.1 -O3 -c -module modules Definitions.f90 -o Definitions.o
ifort -axsse4.1 -O3 -c -module modules Read_input_data.f90 -o Read_input_data.o
ifort -axsse4.1 -O3 -c -module modules Main.f90 -o Main.o
ifort -axsse4.1 -O3 -c -module modules InitializeCheckerboard.f90 -o InitializeCheckerboard.o
ifort -axsse4.1 -O3 -c -module modules MCsweep_Checkerboard.f90 -o MCsweep_Checkerboard.o
ifort -axsse4.1 -O3 -c -module modules Shift_Cells.f90 -o Shift_Cells.o
ifort -axsse4.1 -O3 -c -module modules CB_system_energy_gpu.f90 -o CB_system_energy_gpu.o
ifort -axsse4.1 -O3 -c -module modules energy.f90 -o energy.o
ifort -axsse4.1 -O3 -c -module modules Volume_move.f90 -o Volume_move.o
ifort -axsse4.1 -O3 -c -module modules UpdateCheckerboard.f90 -o UpdateCheckerboard.o

nvcc -ccbin icpc -c Subsweep_Energy_CUDA.cu -o Subsweep_Energy_CUDA.o

ifort  -O3 -mkl -L/usr/local/cuda/lib64 -o bin/MC_GPU_NpT_github.exe Definitions.o Main.o Read_input_data.o InitializeCheckerboard.o MCsweep_Checkerboard.o Shift_Cells.o Subsweep_Energy_CUDA.o CB_system_energy_gpu.o energy.o Volume_move.o UpdateCheckerboard.o -lcuda -lcudart -lstdc++ -lcurand


Input files:

The code needs three input files:
-input.d -> it contains details of the simulations (MC cycles, temperarture, etc..) and of the model potential (number of particle types, number and position of the patches, etc.)
-coords-0.1.dat -> it contains the initial configuration (first line: Number of atoms (N), second-fourth line: vectors defining the simulation box, followed by N lines with four columns: particle label, and Cartesian coordinates of the particle position, followed with other N lines with four entries containing the quaternion that represents the orientation of a given particle).
-the code asks by screen the GPU id number to use
