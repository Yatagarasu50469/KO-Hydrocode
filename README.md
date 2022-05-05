# GENERAL INFORMATION

	**Program Name:**	KO_FiniteDifference
	**Author:**		David Helminiak
	**Version:**		0.2.0
	**License:**		GNU General Public License v3.0
	**Versioning:**		
				0.1.0   Single-threaded re-implementation of finite difference method
				0.2.0   Integration of CPU multiprocessing acceleration through OpenMP

# DESCRIPTION
This hydrocode was ultimately derived from Wilkin's 'Computer Simulation of Dynamic Phenomena.' More specifically, it is a re-implementation of an existing single-threaded version produced in MATLAB by Nathaniel S. Helminiak, itself built on a FORTRAN implementation made available by John P. Borg of the Marquette Unviersity Shock Physics Laboratory. 

# COMPILATION

    mpicc DH_KO_FiniteDifference\ -\ singleThreadBackup.c -o KO_FiniteDifference -lm

# RUNTIME
Change number of threads to parallelize operations across with:
	export OMP_NUM_THREADS=8
	
Then run with:
    ./KO_FiniteDifference

