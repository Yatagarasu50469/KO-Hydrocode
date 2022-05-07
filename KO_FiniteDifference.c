//=================================================================================================================
//Program Name:		KO_FiniteDifference
//Author:			David Helminiak
//Version:			0.2.0
//Date Created:		3 April 2022
//Date Modified: 	7 May 2022
//
//License:			GNU General Public License v3.0
//Purpose:			Simulation of materials for one space dimension and time
//Compilation:		mpicc KO_FiniteDifference.c -o KO_FiniteDifference -lm -fopenmp
//Runtime:			./KO_FiniteDifference
//					Change number of threads with export OMP_NUM_THREADS=8
//Notes:			This hydrocode was ultimately derived from Wilkin's 'Computer Simulation of Dynamic Phenomena.' 
//					More specifically, it is a re-implementation of an existing single-threaded version produced
//                  in MATLAB by Nathaniel S. Helminiak, itself built on a FORTRAN implementation made available
//					by John P. Borg of the Marquette Unviersity Shock Physics Laboratory. 
//
//Versioning:		0.1.0   Single-threaded re-implementation of finite difference method
//					0.2.0	Integration of CPU multiprocessing acceleration through OpenMP
//
//=================================================================================================================

//=================================================================================================================
//External Libraries and Global Definitions
//=================================================================================================================

#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

//=================================================================================================================

//=================================================================================================================
//Function Definitions
//=================================================================================================================



//=================================================================================================================


//=================================================================================================================
//Main Program
//=================================================================================================================

int main(int argc, char **argv) {
	
	//Start a clock to store computation time
	struct timeval startTime;
	struct timeval endTime;
	gettimeofday(&startTime, 0);
	
	//=================================================================================================================
	//CONFIGURATION
	//=================================================================================================================
	
	//Simulation parameters
	//=================================================================================================================
	
	//Time to stop simulation (us)
	int tstop=10;
	
	//Geometry type; plane: 1; cylindircal: 2 (untested); sphere: 3 (untested)
	int geometry = 1;
	
	//Artificial viscosity constant 1; Default: 2
	double artViscConstant1 = 2;
	
	//Artificial viscosity constant 2; Default: 1
	double artViscConstant2 = 1;

	//Pressure inside void positions
	double voidPressure = 0;
	
	//Timestep
	double deltat = 3e-5;

	//Adaptive mesh refinement
	double adaptiveMeshRefinement = 0.0;
	
	//Definition of simulation materials [Number of materials][Number of material properties] in the simulation
		//	0, 		1, 			2, 			3, 			4, 				5, 			6, 			7, 			8, 			9, 			10, 	11, 	12, 		13, 			14, 		15, 		16, 		17, 			18, 	19, 	20, 	21    
		//	ieos, 	Nodes, 		Length, 	xstart, 	P_0,			U_0, 		Rho_0, 		E_0, 		Rho_00, 	Co, 		s,		s2, 	gamma,		Yield (Y), 		mu (G),		pfrac, 		Cv, 		a_p, 			a_e, 	p_e, 	p_s, 	rho_0  
		//						cm			cm			Mabar			cm/us		g/cc					g/cc 		cm/us									Mbar			Mbar		Mbar		Mbar-cc/K-g
	
	
	//CONFIG-A
	
	
	//1x No-spall condition simulation (0: ZEROS; 1: Tungsten; 2: Aluminium; 3: Aluminium; 4: PMMA) (680 material nodes)
	double matProp[5][22] = {
		{	0,    	0,    		0, 	  		0, 			0, 				0,    		0,  		0,	 	    0, 	 	  	0,  	   	0, 	   	0,  	0, 			0, 	  			0, 	 		0, 			0,   		0,   			0, 	  	0,   	0,   	0		},
		{	1, 		320,		8.0,		-8.01,		0.0*1e-6, 		1650*1e-4,	19.20, 		0.000, 		19.20, 		0.404, 		1.23, 	0.00, 	0.00,		1.5*1e-99, 		.35, 		1e9, 		1e-9, 		0.00, 			0.00, 	0.00, 	0.00, 	0.00	},
		{	1,  	20,		   	0.5,     	0.00,   	0.0*1e-6, 		1650*1e-4,  2.71,  		0.000,    	2.71,    	0.535,   	1.34,   0.00,   0.00,   	1.5*1e-99,  	.33,  		1e9,    	1e-9,   	0.00,   		0.00,   0.00,   0.00,   0.00	},
		{	1, 		20, 		0.5, 		0.50, 		0.0*1e-6, 		0.00000, 	2.71, 		0.000, 		2.71, 		0.535, 		1.34, 	0.00, 	0.00, 		1.5*1e-99, 		.33, 		1e9, 		1e-9, 		0.00, 			0.00, 	0.00, 	0.00, 	0.00	},
		{	1, 		320, 		8.0, 		1.00, 		0.0*1e-6, 		0.00000, 	1.19, 		0.000, 		1.19, 		0.257, 		1.54, 	0.00, 	0.00, 		1.5*1e-99, 		.33, 		1e9, 		1e-9, 		0.00, 			0.00, 	0.00, 	0.00, 	0.00	}
	};
	
	
	/*
	//2x No-spall condition simulation (0: ZEROS; 1: Tungsten; 2: Aluminium; 3: Aluminium; 4: PMMA) (1360 material nodes)
	double matProp[5][22] = {
		{	0,    	0,    		0, 	  		0, 			0, 				0,    		0,  		0,	 	    0, 	 	  	0,  	   	0, 	   	0,  	0, 			0, 	  			0, 	 		0, 			0,   		0,   			0, 	  	0,   	0,   	0		},
		{	1, 		640,		16.0,		-16.01,		0.0*1e-6, 		1650*1e-4,	19.20, 		0.000, 		19.20, 		0.404, 		1.23, 	0.00, 	0.00,		1.5*1e-99, 		.35, 		1e9, 		1e-9, 		0.00, 			0.00, 	0.00, 	0.00, 	0.00	},
		{	1,  	40,		   	1.0,     	0.00,   	0.0*1e-6, 		1650*1e-4,  2.71,  		0.000,    	2.71,    	0.535,   	1.34,   0.00,   0.00,   	1.5*1e-99,  	.33,  		1e9,    	1e-9,   	0.00,   		0.00,   0.00,   0.00,   0.00	},
		{	1, 		40, 		1.0, 		1.00, 		0.0*1e-6, 		0.00000, 	2.71, 		0.000, 		2.71, 		0.535, 		1.34, 	0.00, 	0.00, 		1.5*1e-99, 		.33, 		1e9, 		1e-9, 		0.00, 			0.00, 	0.00, 	0.00, 	0.00	},
		{	1, 		640, 		16.0, 		2.00, 		0.0*1e-6, 		0.00000, 	1.19, 		0.000, 		1.19, 		0.257, 		1.54, 	0.00, 	0.00, 		1.5*1e-99, 		.33, 		1e9, 		1e-9, 		0.00, 			0.00, 	0.00, 	0.00, 	0.00	}
	};
	*/
	
	/*
	//3x No-spall condition simulation (0: ZEROS; 1: Tungsten; 2: Aluminium; 3: Aluminium; 4: PMMA) (2040 material nodes)
	double matProp[5][22] = {
		{	0,    	0,    		0, 	  		0, 			0, 				0,    		0,  		0,	 	    0, 	 	  	0,  	   	0, 	   	0,  	0, 			0, 	  			0, 	 		0, 			0,   		0,   			0, 	  	0,   	0,   	0		},
		{	1, 		960,		24.0,		-24.01,		0.0*1e-6, 		1650*1e-4,	19.20, 		0.000, 		19.20, 		0.404, 		1.23, 	0.00, 	0.00,		1.5*1e-99, 		.35, 		1e9, 		1e-9, 		0.00, 			0.00, 	0.00, 	0.00, 	0.00	},
		{	1,  	60,		   	1.5,     	0.00,   	0.0*1e-6, 		1650*1e-4,  2.71,  		0.000,    	2.71,    	0.535,   	1.34,   0.00,   0.00,   	1.5*1e-99,  	.33,  		1e9,    	1e-9,   	0.00,   		0.00,   0.00,   0.00,   0.00	},
		{	1, 		60, 		1.5, 		1.50, 		0.0*1e-6, 		0.00000, 	2.71, 		0.000, 		2.71, 		0.535, 		1.34, 	0.00, 	0.00, 		1.5*1e-99, 		.33, 		1e9, 		1e-9, 		0.00, 			0.00, 	0.00, 	0.00, 	0.00	},
		{	1, 		960, 		24.0, 		3.00, 		0.0*1e-6, 		0.00000, 	1.19, 		0.000, 		1.19, 		0.257, 		1.54, 	0.00, 	0.00, 		1.5*1e-99, 		.33, 		1e9, 		1e-9, 		0.00, 			0.00, 	0.00, 	0.00, 	0.00	}
	};
	*/
	
	
	//CONFIG-B
	/*
	//1x Spall condition simulation (0: ZEROS; 1: Tungsten; 2: Aluminium) (680 material nodes)
	double matProp[3][22] = {
		{	0,    	0,    		0, 	  		0, 			0, 				0,    		0,  		0,	 	    0, 	 	  	0,  	   	0, 	   	0,  	0, 			0, 	  			0, 	 		0, 			0,   		0,   			0, 	  	0,   	0,   	0		},
		{	1, 		640,		16.0,		-16.01,		0.0*1e-6, 		1650*1e-4,	19.20, 		0.000, 		19.20, 		0.404, 		1.23, 	0.00, 	0.00,		1.5*1e-99, 		.35, 		1e9, 		1e-9, 		0.00, 			0.00, 	0.00, 	0.00, 	0.00	},
		{	1,  	40,		   	1.0,     	0.00,   	0.0*1e-6, 		0.00000, 	2.71,  		0.000,    	2.71,    	0.535,   	1.34,   0.00,   0.00,   	1.5*1e-99,  	.33,  		1e9,    	1e-9,   	0.00,   		0.00,   0.00,   0.00,   0.00	},
	};
	*/
	/*
	//2x Spall condition simulation (0: ZEROS; 1: Tungsten; 2: Aluminium) (1360 material nodes)
	double matProp[3][22] = {
		{	0,    	0,    		0, 	  		0, 			0, 				0,    		0,  		0,	 	    0, 	 	  	0,  	   	0, 	   	0,  	0, 			0, 	  			0, 	 		0, 			0,   		0,   			0, 	  	0,   	0,   	0		},
		{	1, 		1280,		32.0,		-32.01,		0.0*1e-6, 		1650*1e-4,	19.20, 		0.000, 		19.20, 		0.404, 		1.23, 	0.00, 	0.00,		1.5*1e-99, 		.35, 		1e9, 		1e-9, 		0.00, 			0.00, 	0.00, 	0.00, 	0.00	},
		{	1,  	80,		   	2.0,     	0.00,   	0.0*1e-6, 		0.00000, 	2.71,  		0.000,    	2.71,    	0.535,   	1.34,   0.00,   0.00,   	1.5*1e-99,  	.33,  		1e9,    	1e-9,   	0.00,   		0.00,   0.00,   0.00,   0.00	},
	};
	*/
	/*
	//3x Spall condition simulation (0: ZEROS; 1: Tungsten; 2: Aluminium) (2040 material nodes)
	double matProp[3][22] = {
		{	0,    	0,    		0, 	  		0, 			0, 				0,    		0,  		0,	 	    0, 	 	  	0,  	   	0, 	   	0,  	0, 			0, 	  			0, 	 		0, 			0,   		0,   			0, 	  	0,   	0,   	0		},
		{	1, 		1920,		48.0,		-48.01,		0.0*1e-6, 		1650*1e-4,	19.20, 		0.000, 		19.20, 		0.404, 		1.23, 	0.00, 	0.00,		1.5*1e-99, 		.35, 		1e9, 		1e-9, 		0.00, 			0.00, 	0.00, 	0.00, 	0.00	},
		{	1,  	120,		3.0,     	0.00,   	0.0*1e-6, 		0.00000, 	2.71,  		0.000,    	2.71,    	0.535,   	1.34,   0.00,   0.00,   	1.5*1e-99,  	.33,  		1e9,    	1e-9,   	0.00,   		0.00,   0.00,   0.00,   0.00	},
	};
	*/
	
	//Set boundary conditions
	//		0		1			2		3					4			5
	//		iend	P_0			U_0		Rho_0				E_0			V_0
	double boundary[2][6] = {
		{	-1, 	0.0*1e-1, 	0.0, 	matProp[0][8], 		0.000, 		1.0 	},
		{	1,		0.0*1e-1, 	0.0, 	matProp[1][8],		0.000, 		1.0 	}
	};
	
	
	//Operational parameters
	//=================================================================================================================
	
	//Number of steps between printing
	int printIterSkip = 1;
	
	//Step to begin writing data to disk
	int fileWrite = 0;
	
	//Time between writing files
	double fileTimeSkip = 0.01;
	
	//Enable/Disable debug mode
	bool debug = false; 
	
	//If debug mode is enabled, specify the node to examine to print information about
	int debugNode = 2;
	
	//=================================================================================================================
	//SETUP
	//=================================================================================================================
	
	//Static program parameters
	//=================================================================================================================
	
	//Initialize adjustable functional timestep and a duplicate as a alternative when no node contact occurs
	double delt = deltat;
	double delt2 = deltat;
	
	//Current timestep
	int timestep = 0;
	
	//Indicate if contact occurs (1)
	bool contactFlag = false;
	
	//Indicate if timestep is adjusted
	bool timestepAdjustFlag = false;

	//Number of defined materials (-1 for zero material)
	int numMat = sizeof(matProp)/sizeof(matProp[0]);
	
	//Number of allowed material nodes; at least double the number of starting nodes for fracture 
	int totalNodes = 0;
	#pragma omp parallel for reduction(+:totalNodes)
	for (int i=0; i<numMat; i++) {
		totalNodes += matProp[i][1];
	}
	totalNodes = totalNodes*2;
	
	//Initialize variables
	//=================================================================================================================
	
	//Pressure of each node across applicable time steps
	double P[4][totalNodes];
	
	//Mass of each node
	double m[totalNodes];
	
	//Energy of each node across applicable time steps
	double E[4][totalNodes];
	
	//Volume of each node across applicable time steps
	double V[4][totalNodes];
	
	//Density of each node
	double rho[totalNodes];
	
	//Initial boundary condition?; 9 specifies unused nodes
	int ibc[totalNodes];
	
	//Time references (Past: 0, Current: 1, Next: 2, 3)
	double t[4];
	
	//
	double Y[totalNodes];
	
	//
	double pfrac[totalNodes];
	
	//Velocities of each node across applicable time steps
	double U[4][totalNodes];
	double U2[totalNodes];
	
	//Positions of each node
	double r[4][totalNodes];
	
	//
	double phi[4][totalNodes];
	//
	
	double sigmar[4][totalNodes];
	double sigmao[4][totalNodes];
	
	//
	double beta[4][totalNodes];
	
	//
	double q[4][totalNodes];
	
	//
	double s1[4][totalNodes];
	double s2[4][totalNodes];
	double s3[4][totalNodes];
	
	//
	double epsilon1[4][totalNodes];
	double epsilon2[4][totalNodes];
	
	//
	double K[4][totalNodes];
	
	//Temperature of each node across applicable time steps (K)
	double T[4][totalNodes];
	
	//Entropy of each node across applicable time steps (mbar-cc/K/g)
	double entropy[4][totalNodes];

	//Initial equations of state for each node
	int ieos[2][totalNodes];
	
	//
	int icompact[4][totalNodes];
	
	//Holding variable for velocity
	double uSave;
	
	//Holding variables for temporary resolution of inside/outside boundary conditions
	double sigmar_in, sigmao_in, sigmar_out, sigmao_out;
	
	//Contact times for resolution of node-to-node contact
	double contactTime, initialContactTime;
	
	//Inside and outside boundary conditon variables for resolution of node-to-node contact
	double phiv, phij, betav, betaj;
	
	//Other unlabeled variables that are called or referenced...
	double deltaZ, Vdot, dt_min, delt_temp;
	double qbar;
	double b;
	double aa, bb, xb;
	double xx, xa;
	double qtotal, mvtotal, ketotal, ietotal, etotal, ke3total;
	double bs1, bs2, bs3, bs4, bs5, bs6, bs7, bs8, bs9, bs10;
	double rho_local, stemp, ctemp, gtemp, v0, v00, vv;
	double k1, k2, k3, gamma0;
	double Us, up, PH, EH, TH, strain, P0, E0, T0;
	double En2j1, diffE;
	
	//Generic variable to hold a total value
	double total;
	
	//Generic variable to use as a loop counter as needed
	int counter;
	
	//Declare index variable for material nodes
	int iPoint = 0;
	
	//A material's length/# of material nodes (half those assigned)
	double deltar;
	
	//In the event of node-to-node contact, the minimum found contact time determined
	double minContactTime = -1;
	
	//In the event of node-to-node contact, the node found corresponding to the minimum found contact time
	int contactNode;
	
	//In the event that an error occurs inside of a parallelized section, trigger program termination
	bool errorCode = false;
	
	//Specify initial values, confirming an even number of material nodes
	#pragma omp parallel for
	for (int j=0; j<totalNodes; j++) {
		m[j] = 0;
		rho[j] = 0;
		ibc[j] = 9; //Specify all nodes as unused to start
		Y[j] = 0;
		U2[j] = 0;
		pfrac[j] = 0;
		ieos[0][j] = 0;
		ieos[1][j] = 0;
		for (int n=0; n<4; n++) {
			P[n][j] = 0;
			E[n][j] = 0;
			V[n][j] = 0;
			U[n][j] = 0;
			r[n][j] = 0;
			phi[n][j] = 0;
			sigmar[n][j] = 0;
			sigmao[n][j] = 0;
			beta[n][j] = 0;
			q[n][j] = 0;
			s1[n][j] = 0;
			s2[n][j] = 0;
			s3[n][j] = 0;
			epsilon1[n][j] = 0;
			epsilon2[n][j] = 0;
			K[n][j] = 0;
			T[n][j] = 0;
			entropy[n][j] = 0;
			icompact[n][j] = 0;
		}
	}
	
	#pragma omp parallel for
	for (int n=0; n<4; n++) {
		t[n] = deltat * (double) n;
	}
	
	//Check for an even number of material nodes
	#pragma omp parallel for
	for (int imat=1; imat<numMat; imat++) {
		ieos[0][totalNodes-(imat-1)] = imat;
		ieos[1][totalNodes-(imat-1)] = matProp[imat][0];
		if ((int)matProp[imat][1]%2 != 0) {
			printf("Error - Number of material nodes must be even, please check parameters for material: %d\n", imat);
			errorCode = true;
		}
	}
	
	//Set initial boundaries at first and last nodes; (-/+) Free end: 1, Contact check: 2 and 3, Fixed end: 4
	ibc[0] = boundary[0][0];
	ibc[totalNodes] = boundary[1][0];
	
	//Discretization
	//=================================================================================================================

	//Obtain initial positions from the material properties table
	r[0][1] = matProp[0][3];
	r[1][1] = matProp[0][3];
	
	//Discretize nodes
	for (int imat=1; imat<numMat; imat++) {
		
		if (imat == 1) {
			iPoint = 0;
		}
		//If there is no gap between materials
		else if (fabs(r[0][iPoint+(int)matProp[imat-1][1]]-matProp[imat][3]) < 1e-5) {
			
			//Assign initial positions and boundary conditions and move node pointer to the next material position
			r[0][iPoint+(int)matProp[imat-1][1]] = matProp[imat][3];
			r[1][iPoint+(int)matProp[imat-1][1]] = matProp[imat][3];
			ibc[iPoint+(int)matProp[imat-1][1]] = 0;
			iPoint = iPoint+matProp[imat-1][1];
		}
		//If there is an initial gap between materials
		else if (r[0][iPoint+(int)matProp[imat-1][1]] < matProp[imat][3]) {
			
			//Advance considered node; skipping fracture node
			iPoint = iPoint+matProp[imat-1][1]+2;
			
			//Assign boundary conditions for fracture and material nodes
			ibc[iPoint] = -2;
			ibc[iPoint-2] = 2;
			
			//Set initial equation of state variable
			ieos[1][iPoint-1] = 0;
		}
		else {
		   printf("Error - Input geometry issue encountered during discretization of material %d\n", imat);
		   return -1;
		   
		}
		
		//Nodes are defined in pairs, where the second is used in case of fracture; set parameters for the first node pair
		r[0][iPoint] = matProp[imat][3];
		r[1][iPoint] = matProp[imat][3];
		U[0][iPoint] = matProp[imat][5];
		U[1][iPoint] = matProp[imat][5];
		deltar = matProp[imat][2]/(matProp[imat][1]/2);
		
		//Define each remaining node pair
		#pragma omp parallel for private(counter)
		for (int j=iPoint+2; j<=iPoint+(int)matProp[imat][1]; j+=2) {
			counter = (j/2)-(iPoint/2);
			ibc[j] = 0;
			ibc[j-1] = 0;
			r[0][j] = matProp[imat][3]+(counter*deltar);
			r[1][j] = matProp[imat][3]+(counter*deltar);
			r[0][j-1] = 0.5*(r[0][j]+matProp[imat][3]+((counter-1)*deltar));
			r[1][j-1] = 0.5*(r[1][j]+matProp[imat][3]+((counter-1)*deltar));
			U[0][j] = matProp[imat][5];
			U[1][j] = matProp[imat][5];
			P[0][j-1] = matProp[imat][4];
			P[1][j-1] = matProp[imat][4];
			E[0][j-1] = matProp[imat][7];
			E[1][j-1] = matProp[imat][7];
			ieos[0][j-1] = ieos[0][totalNodes-(imat-1)];
			ieos[1][j-1] = ieos[1][totalNodes-(imat-1)];
			rho[j-1] = matProp[imat][6];  //rho is always rho0; DH, so why not replace rho entirely with a matProp call...
			Y[j-1] = matProp[imat][13];			
			pfrac[j-2] = matProp[imat][15];
			pfrac[-1] = matProp[imat][15];
			pfrac[j] = matProp[imat][15];
		}
		
		//Assign values to the last cell center
		ieos[0][totalNodes-(imat-1)] = 0;
		ieos[1][totalNodes-(imat-1)] = 0;
		ibc[(int)matProp[imat][1]+iPoint] = 2;
		ibc[(int)matProp[imat][1]+iPoint+1] = 9;
	}
	
	//Assign values to the last node
	ibc[(int)matProp[numMat-1][1]+iPoint] = ibc[totalNodes];
	ibc[totalNodes]=9;

	//Assign mass, volume, energy, temperature, and entropy to non-fracture nodes
	#pragma omp parallel for
	for (int j=0; j<=totalNodes-2; j+=2) {
		if (ibc[j+1] == 0) {
			m[j+1] = rho[j+1]*((pow(r[0][j+2],geometry)-pow(r[0][j],geometry))/geometry);
			V[0][j+1] = rho[j+1]*((pow(r[0][j+2],geometry)-pow(r[0][j],geometry))/geometry)/m[j+1];
			V[1][j+1] = rho[j+1]*((pow(r[0][j+2],geometry)-pow(r[0][j],geometry))/geometry)/m[j+1];
			E[0][j+1] = P[0][j+1]/(V[0][j+1]*matProp[ieos[0][j+1]][12]-1);
			E[1][j+1] = P[1][j+1]/(V[1][j+1]*matProp[ieos[0][j+1]][12]-1);
			T[0][j+1] = E[0][j+1]/(rho[j+1]*matProp[ieos[0][j+1]][16]);
			T[1][j+1] = E[1][j+1]/(rho[j+1]*matProp[ieos[0][j+1]][16]);
			entropy[0][j+1] = 6.8*1e-5;
			entropy[1][j+1] = 6.8*1e-5;
			
			//#pragma omp critical
			//{
			//printf("Node: %d\n", j);
			//printf("Pressure (Mbar): %2.0f\n", P[0][j+1]);
			//printf("Energy (MBar-cc/cc): %2.0f\n", E[0][j+1]);
			//printf("Volume (cc): %2.0f\n", V[0][j+1]);
			//printf("Temperature (K): %2.0f\n", T[0][j+1]);
			//printf("Density (g/cc): %2.0f\n", rho[j+1]);
			//printf("Mass (g): %2.0f\n", m[j+1]);
			//printf("Entropy: %2.0f\n", entropy[0][j+1]);
			//printf("\n");
			//}
		}
	}
	
	//Perform simulation until the current time reaches the stopping criteria
	while (t[1] <= tstop) {
		
		//Reset the contact and timestep adjustment flags
		contactFlag = false;
		timestepAdjustFlag = false;
		
		//Advance the current time step
		timestep += 1;
		
		//Contact Checks and Time Advancement
		//=================================================================================================================
		
		//Adjust future node positions according to velocities
		#pragma omp parallel for
		for (int j=0; j<=totalNodes; j+=2) {
			if (ibc[j] != 9) {
				r[3][j] = r[1][j]+U[0][j]*deltat;
				r[2][j] = (r[1][j]+r[3][j])/2.0;
			}
		}
		
		//Perform a contact check for non-void nodes
		#pragma omp parallel for private(sigmar_in, sigmao_in, sigmar_out, sigmao_out, phiv, betav, phij, betaj, aa, bb, counter, contactTime, initialContactTime)
		for (int j=2; j<=totalNodes; j+=2) {
			
			//If the current mode is not void, its next position intersects with the previous node, and the fracture node is void
			if ((ibc[j] != 9) && (r[3][j] <= r[3][j-2]) && (ibc[j-1] == 9)) {
				
				//Resolve inside boundary condition
				sigmar_in = -(P[1][j-3]+q[0][j-3])+s1[1][j-3];
				sigmao_in = -(P[1][j-3]+q[0][j-3])+s2[1][j-3];
				phiv = (rho[j-3]*((r[1][j-2]-r[1][j-4])/V[1][j-3]))/2.0;
				betav = (sigmar_in-sigmao_in)*V[1][j-1]/(r[1][j-3]*rho[j-3]);
				
				//Resolve outside boundary condition
				sigmar_out = (-(P[1][j+1]+q[0][j+1])+s1[1][j+1]);
				sigmao_out = (-(P[1][j+1]+q[0][j+1])+s2[1][j+1]);
				phij = (rho[j+1]*((r[1][j+2]-r[1][j])/V[1][j+1]))/2.0;
				betaj = (sigmar_out-sigmao_out)*V[1][j+1]/(r[1][j+1]*rho[j+1]);
				
				//Compute time until nodes contact
				aa = sigmar_out/phij+sigmar_in/phiv+(betav+betaj)*(geometry-1.0);
				bb = (2.0*(U[0][j]-U[0][j-2])+aa*(t[1]-t[0]));
				contactTime = -(2.0*(r[1][j]-r[1][j-2]))/(bb);
				
				//Adjust timestep until the time between contacting nodes is small enough to ignore, or a maximum number of iterations is reached
				counter = 0;
				initialContactTime = 0.0;
				while ((fabs(contactTime-initialContactTime) > 1e-16) && counter < 10*1e3) {
					counter+=1;
					initialContactTime = contactTime;
					contactTime = contactTime-(aa*contactTime*contactTime+bb*contactTime+2.0*(r[1][j]-r[1][j-2]))/(2.0*aa*contactTime+bb);
				}
				
				#pragma omp critical
				{
					//Indicate occurrence of future contact 
					contactFlag = true;
					
					//If the minimum contact time has not yet been set, or the determined time for this node contact is less than the prior, then update such
					if ((minContactTime == -1) || (contactTime < minContactTime)) {
						minContactTime = contactTime;
						contactNode = j;
						
						//DH, not sure if setting sigmar and sigmao is actually needed here, or if the computation is repeated/overwrites this later on
						//If repeated, then can replace the sigmar/o-in/out holding variables, otherwise the access needs to be protected
						sigmar[1][j-3] = sigmar_in;
						sigmao[1][j-3] = sigmao_in;
						sigmar[1][j+1] = sigmar_out;
						sigmao[1][j+1] = sigmao_out;
					}
				}
			}
		}
		
		//If contact, set timestep/boundary condtion to the minimum contact time/node and resolve in/out boundary condition; either way compute new timesteps
		if (contactFlag) {
			deltat = minContactTime/2.0;
			ibc[contactNode] = -3;
			t[2] = t[1] + minContactTime/2.0;
			t[3] = t[2] + deltat;
		}
		else {
			t[2] = t[1] + deltat/2.0;
			t[3] = t[1] + deltat;
		}
		
		//Compute resulting default (and alternate if no contact) functional time difference(s)
		delt = t[2]-t[0];
		if (!contactFlag) { 
			delt2 = t[1]+(deltat/8)-t[0];
		}
		
		//Ensure Conservation of Momentum
		//=================================================================================================================
		
		//Update sigmar/sigmao for nodes with void boundary conditions
		#pragma omp parallel for
		for (int j=0; j<=totalNodes; j+=2) {
			if ((j < totalNodes) && (ibc[j+1] == 0)) {
				sigmar[1][j+1] = -(P[1][j+1]+q[0][j+1])+s1[1][j+1];
				sigmao[1][j+1] = -(P[1][j+1]+q[0][j+1])+s2[1][j+1];
			}
		}

		//For nodes that are not voids, adjust velocities depending on the initial boundary conditions
		#pragma omp parallel for
		for (int j=0; j<=totalNodes; j+=2) {
			if (ibc[j] != 9) {
				if (ibc[j] == -1) {
					phi[1][j] = (rho[j+1]*((r[1][j+2]-r[1][j])/V[1][j+1]))/2.0;
					beta[1][j] = (sigmar[1][j+1]-sigmao[1][j+1])*V[1][j+1]/(r[1][j+1]*rho[j+1]);
					if ((delt/phi[1][j])*sigmar[1][j+1]+delt*beta[1][j]*(geometry-1.0) == 0) {
						U[2][j] = U[0][j];
					}
					else {
						U[2][j] = U[0][j]+(delt/phi[1][j])*sigmar[1][j+1]+delt*beta[1][j]*(geometry-1.0);
					}
					if (!contactFlag) {
						if ((delt2/phi[1][j])*sigmar[1][j+1]+delt2*beta[1][j]*(geometry-1.0) == 0) {
							U2[j] = U[0][j];
						}
						else {
							U2[j] = U[0][j]+(delt2/phi[1][j])*sigmar[1][j+1]+delt2*beta[1][j]*(geometry-1.0);
						}
					}
				}
				else if (ibc[j] == 1) {
					phi[1][j] = rho[j-1]*((r[1][j]-r[1][j-2])/V[1][j-1])/2.0;
					beta[1][j] = (sigmar[1][j-1]-sigmao[1][j-1])*V[1][j-1]/(r[1][j-1]*rho[j-1]);
					if ((delt/phi[1][j])*(-boundary[1][2]-sigmar[1][j-1])+delt*beta[1][j]*(geometry-1.0) == 0) {
						U[2][j] = U[0][j];
					}
					else {
						U[2][j] = U[0][j]+(delt/phi[1][j])*(-boundary[1][2]-sigmar[1][j-1])+delt*beta[1][j]*(geometry-1.0);
					}
					if (!contactFlag) {
						if ((delt2/phi[1][j])*(-boundary[1][2]-sigmar[1][j-1])+delt2*beta[1][j]*(geometry-1.0) == 0) {
							U2[j] = U[0][j];
						}
						else {
							U2[j] = U[0][j]+(delt2/phi[1][j])*(boundary[1][2]-sigmar[1][j-1])+delt2*beta[1][j]*(geometry-1.0);
						}
					}
				}
				else if (ibc[j] == -2 || ibc[j] == -3) {
					phi[1][j] = (rho[j+1]*((r[1][j+2]-r[1][j])/V[1][j+1]))/2.0;
					beta[1][j] = (sigmar[1][j+1]-sigmao[1][j+1])*V[1][j+1]/(r[1][j+1]*rho[j+1]);
					U[2][j] = U[0][j]+(delt/phi[1][j])*(sigmar[1][j+1]+voidPressure)+delt*beta[1][j]*(geometry-1.0);
					if (!contactFlag) {
						U2[j] = U[0][j]+(delt2/phi[1][j])*(sigmar[1][j+1]+voidPressure)+delt2*beta[1][j]*(geometry-1.0);
					}
				}
				else if (ibc[j] == 2) {
					phi[1][j] =  rho[j-1]*((r[1][j]-r[1][j-2])/V[1][j-1])/2.0;
					beta[1][j] = (sigmar[1][j-1]-sigmao[1][j-1])*V[1][j-1]/(r[1][j-1]*rho[j-1]);
					U[2][j] = U[0][j]+(delt/phi[1][j])*(voidPressure-sigmar[1][j-1])+delt*beta[1][j]*(geometry-1.0);
					if (!contactFlag) {
						U2[j] = U[0][j]+(delt2/phi[1][j])*(voidPressure-sigmar[1][j-1])+delt2*beta[1][j]*(geometry-1.0);
					}
				}
				else if (ibc[j] == 0) {
					phi[1][j] = (0.5)*(rho[j+1]*((r[1][j+2]-r[1][j])/V[1][j+1])+rho[j-1]*((r[1][j]-r[1][j-2])/V[1][j-1]));
					beta[1][j]=((sigmar[1][j+1]-sigmao[1][j+1])*(V[1][j+1]/rho[j+1])/(0.5*(r[1][j+2]+r[1][j]))+(sigmar[1][j-1]-sigmao[1][j-1])*(V[1][j-1]/rho[j-1])/(0.5*(r[1][j]+r[1][j-2])))/2.0;
					if ((delt/phi[1][j])*(sigmar[1][j+1]-sigmar[1][j-1])+delt*beta[1][j]*(geometry-1.0) == 0) {
						U[2][j] = U[0][j];
					}
					else {
						U[2][j] = U[0][j]+(delt/phi[1][j])*(sigmar[1][j+1]-sigmar[1][j-1])+delt*beta[1][j]*(geometry-1.0);
					}
					if (!contactFlag) {
						if ((delt2/phi[1][j])*(sigmar[1][j+1]-sigmar[1][j-1])+delt2*beta[1][j]*(geometry-1.0) == 0) {
							U2[j] = U[0][j];
						}
						else {
							U2[j] = U[0][j]+(delt2/phi[1][j])*(sigmar[1][j+1]-sigmar[1][j-1])+delt2*beta[1][j]*(geometry-1.0);
						}
					}
				}
				else if (abs(ibc[j]) == 4) {
					U[2][j] = 0.0;
					if (!contactFlag) {
						U2[j] = 0.0;
					}
				}
				else {
					printf("Error - Momentum check failed at node: %d with initial boundary condition: %d\n", j, ibc[j]);
					errorCode = true;
				}
			}
		}
		
		//If contact did not occur then adjust the timestep and compute non-void node velocities, otherwise computations were based on the contact time
		if (!contactFlag) {
			
			//Compute condition for potential timestep adjustment
			counter = 0;
			total = 0.0;
			#pragma omp parallel for reduction(+:total) reduction(+:counter)
			for (int j=0; j<totalNodes; j++) {
				total += U[2][j]-U2[j];
				if (U[2][j] > 0) {
					counter+=1;
				}
			}
			
			//Adjust the true timestep, increasing it if condition is satisfied and otherwise decreasing it and using the values just computed
			if (total/counter < 1e-9) {
				deltat = deltat*1.01;
				timestepAdjustFlag = true;
			}
			else {
				deltat = deltat/4;
				#pragma omp parallel for
				for (int j=0; j<=totalNodes; j+=2) {
					if (ibc[j] != 9) {
						if ((ibc[j] == -1) || (ibc[j] == 1) || (ibc[j] == -2 || ibc[j] == -3) || (ibc[j] == 2) || (ibc[j] == 0)) {
							U[2][j] = U2[j];
						}
					}
				}
			}
			
			//Adjust the functional timestep according to the modification made to the true timestep
			t[2] = t[1] + deltat/2.0;
			t[3] = t[1] + deltat;
			delt = t[2]-t[0];
			
			//If not using the values just computed, then recompute the non-void node velocities based on the adjusted timestep depending on the initial boundary conditions
			if (timestepAdjustFlag) {
				#pragma omp parallel for
				for (int j=0; j<=totalNodes; j+=2) {
					if (ibc[j] != 9) {
						if (ibc[j] == -1) {
							if ((delt/phi[1][j])*sigmar[1][j+1]+delt*beta[1][j]*(geometry-1.0) == 0) {
								U[2][j] = U[0][j];
							}
							else {
								U[2][j] = U[0][j]+(delt/phi[1][j])*sigmar[1][j+1]+delt*beta[1][j]*(geometry-1.0);
							}
						}
						else if (ibc[j] == 1) {
							if ((delt/phi[1][j])*(-boundary[1][2]-sigmar[1][j-1])+delt*beta[1][j]*(geometry-1.0) == 0) {
								U[2][j] = U[0][j];
							}
							else {
								U[2][j] = U[0][j]+(delt/phi[1][j])*(-boundary[1][2]-sigmar[1][j-1])+delt*beta[1][j]*(geometry-1.0);
							}
						}
						else if (ibc[j] == -2 || ibc[j] == -3) {
							U[2][j] = U[0][j]+(delt/phi[1][j])*(sigmar[1][j+1]+voidPressure)+delt*beta[1][j]*(geometry-1.0);
						}
						else if (ibc[j] == 2) {
							U[2][j] = U[0][j]+(delt/phi[1][j])*(voidPressure-sigmar[1][j-1])+delt*beta[1][j]*(geometry-1.0);
						}
						else if (ibc[j] == 0) {
							if ((delt/phi[1][j])*(sigmar[1][j+1]-sigmar[1][j-1])+delt*beta[1][j]*(geometry-1.0) == 0) {
								U[2][j] = U[0][j];
							}
							else {
								U[2][j] = U[0][j]+(delt/phi[1][j])*(sigmar[1][j+1]-sigmar[1][j-1])+delt*beta[1][j]*(geometry-1.0);
							}
						}
						else if (abs(ibc[j]) == 4) {
							U[2][j] = 0.0;
						}
						else {
							printf("Error - Momentum check failed at node: %d with initial boundary condition: %d\n", j, ibc[j]);
							errorCode = true;
						}
					}
				}
			}
		}
		
		//Update positions
		//=================================================================================================================
		
		//Move positions of non-void nodes according to velocities/timestep
		#pragma omp parallel for
		for (int j=0; j<totalNodes; j+=2) {
			if (ibc[j] != 9) {
				r[3][j] = r[1][j]+U[2][j]*deltat;
				r[2][j] = (r[1][j]+r[3][j])/2.0;
			}
		}
		
		//If there was contact, connect contacting nodes and shift all other nodes accordingly
		if (contactFlag) {
			
			//Reset the contact node boundary condition
			ibc[contactNode] = 0;
			
			//Output the condition and action to the console
			//printf("Joining Nodes %d and %d\n",contactNode-2,contactNode);
			//printf("Before Joining\n");
			//printf("\tLEFT:\t\tr[1][contactNode-4]=%2.4f;\tU[0][contactNode-4]=%2.4f\n", r[1][contactNode-4], U[0][contactNode-4]);
			//printf("\tVOID:\t\tr[1][contactNode-2]=%2.4f;\tU[1][contactNode-2]=%2.4f\n", r[1][contactNode-2], U[0][contactNode-2]); 
			//printf("\tJ:\t\tr[1][contactNode]=%2.4f;\tU[0][contactNode]=%2.4f\n", r[1][contactNode], U[0][contactNode]);
			//printf("\tRIGHT:\t\tr[1][contactNode+2]=%2.4f;\tU[0][contactNode+2]=%2.4f\n", r[1][contactNode+2], U[0][contactNode+2]);
			
			//Store current node information
			uSave = (m[contactNode+1]*U[2][contactNode]+m[contactNode-3]*U[2][contactNode-2])/(m[contactNode-3]+m[contactNode+1]);
			
			//Shift node information
			//DEV: Contains loop dependencies and cannot be parallelized
			for (int j=contactNode-2; j<totalNodes-2; j++) {
				ibc[j] = ibc[j+2];
				m[j] = m[j+2];
				pfrac[j] = pfrac[j+2];
				ibc[j+1] = ibc[j+3];
				ieos[0][j+1] = ieos[0][j+3];
				ieos[1][j+1] = ieos[1][j+3];
				rho[j+1] = rho[j+3];
				Y[j+1] = Y[j+3];
				pfrac[j+1] = pfrac[j+3];
				for (int nz=1; nz<=3; nz++) {
					r[nz][j] = r[nz][j+2];
					U[nz][j] = U[nz][j+2];
					phi[nz][j] = phi[nz][j+2];
					beta[nz][j] = beta[nz][j+2];
					sigmar[nz][j+1] = sigmar[nz][j+3];
					sigmao[nz][j+1] = sigmao[nz][j+3];
					V[nz][j+1] = V[nz][j+3];
					s1[nz][j+1] = s1[nz][j+3];
					s2[nz][j+1] = s2[nz][j+3];
					s3[nz][j+1] = s3[nz][j+3];
					E[nz][j+1] = E[nz][j+3];
				}
			}
			
			//Restore node information after shift
			U[2][contactNode-2] = uSave;
			pfrac[contactNode-2] = 1e-2;
			
			//Output the condition and action to the console
			//printf("After Joining:\n");
			//printf("\tLEFT:\t\tr[3][contactNode-4]=%2.4f\tU[2][contactNode-4]=%2.4f\n", r[3][contactNode-4], U[2][contactNode-4]);
			//printf("\tVOID:\t\tr[3][contactNode-2]=%2.4f\tU[2][contactNode-2]=%2.4f\n", r[3][contactNode-2], U[2][contactNode-2]);
			//printf("\tJ:\t\tr[3][contactNode]=%2.4f\tU[2][contactNode]=%2.4f\n", r[3][contactNode], U[2][contactNode]);
			//printf("\tRIGHT:\t\tr[3][contactNode+2]=%2.4f\tU[2][contactNode+2]=%2.4f\n", r[3][contactNode+2], U[2][contactNode+2]);
			//printf("\n");
			
			//If contact check fails, output the condition to the console
			if ((ibc[contactNode] != 9) && (r[3][contactNode] < r[3][contactNode-2])) {
				printf("Error - Contact Check Failed; try lowering the timestep\n");
				printf("\tLEFT:\t\tr[3][contactNode-4]=%2.4f\tU[2][contactNode-4]=%2.4f\n", r[3][contactNode-4], U[2][contactNode-4]);
				printf("\tVOID:\t\tr[3][contactNode-2]=%2.4f\tU[2][contactNode-2]=%2.4f\n", r[3][contactNode-2], U[2][contactNode-2]);
				printf("\tJ:\t\tr[3][contactNode]=%2.4f\tU[2][contactNode]=%2.4f\n", r[3][contactNode], U[2][contactNode]);
				printf("\tRIGHT:\t\tr[3][contactNode+2]=%2.4f\tU[2][contactNode+2]=%2.4f\n", r[3][contactNode+2], U[2][contactNode+2]);
				printf("\n");
				errorCode = true;
			}
		}
		
		//Validate boundary conditions are satisfied after the contact check 
		#pragma omp parallel for
		for (int j=2; j<=totalNodes-2; j+=2) {
			if (ibc[j] != 9 && r[3][j] < r[3][j-2]) {
				printf("Error - Contact Check Failed for node: %d, try lowering the timestep\n",j);
				errorCode = true;
			}
		}
		
		//Relative volume, velocity strains, stresses, yield condition, artificial viscosity, energy, pressure, temperature, entropy, and pressure/entropy convergence
		#pragma omp parallel for private(adaptiveMeshRefinement, gamma0, k1, k2, k3, qbar, deltaZ, strain, xa, xb, xx, vv, v0, v00, stemp, ctemp, gtemp, P0, E0, T0, rho_local, up, Us, PH, EH, TH, En2j1, diffE)
		for (int j=0; j<=totalNodes-2; j+=2) {
			if (ibc[j+1] == 0) {
				
				//Relative Volume
				//=================================================================================================================
				r[2][j+1] = (r[2][j]+r[2][j+2])/2.0;
				r[3][j+1] = (r[3][j]+r[3][j+2])/2.0;
				V[3][j+1] = rho[j+1]*((pow(r[3][j+2],geometry)-pow(r[3][j],geometry))/geometry)/m[j+1];
				V[2][j+1] = rho[j+1]*((pow(r[2][j+2],geometry)-pow(r[2][j],geometry))/geometry)/m[j+1];
				if ((V[3][j+1] == 0) || (V[2][j+1] == 0)) {
					printf("Error - Zero volume: %d %2.0f %2.0f %2.0f\n", j+1, r[3][j+2], r[3][j], V[3][j+1]);
					errorCode = true;
				}
				
				//Veclocity Strains; 1st order accurate			
				//=================================================================================================================
				epsilon1[2][j+1] = (U[2][j+2]-U[2][j])/(r[2][j+2]-r[2][j]);
				if (geometry == 1) {
					epsilon2[2][j+1] = 0.0;
				}
				else {
					epsilon2[2][j+1] = (U[2][j+2]+U[2][j])/(r[2][j+2]+r[2][j]);
				}
			
				//Stresses
				//=================================================================================================================
				if (ieos[1][j+1] != 3) {
					s1[3][j+1] = s1[1][j+1] + 2.0*matProp[ieos[0][j+1]][14]*(epsilon1[2][j+1]*deltat-(V[3][j+1]-V[1][j+1])/(3.0*V[2][j+1]));
					s2[3][j+1] = s2[1][j+1] + 2.0*matProp[ieos[0][j+1]][14]*(epsilon2[2][j+1]*deltat-(V[3][j+1]-V[1][j+1])/(3.0*V[2][j+1]));
					s3[3][j+1] = -(s1[3][j+1]+s2[3][j+1]);
					s1[2][j+1] = (s1[3][j+1]+s1[1][j+1])/2.0;
					s2[2][j+1] = (s2[3][j+1]+s2[1][j+1])/2.0;
					s3[2][j+1] = (s3[3][j+1]+s3[1][j+1])/2.0;
				}
				else {
					s1[3][j+1] = 4.0*1.8e-10*epsilon1[2][j+1]/3.0-1.387e-6*(V[3][j+1]-V[1][j+1])/V[2][j+1];
					s1[2][j+1] = (s1[3][j+1]+s1[1][j+1])/2.0;
				}
			
				//Von Mises Yield Condition (Calculate the deviatoric strain at 2 and compare it to the yield strength)
				//=================================================================================================================
				if (ieos[1][j+1] != 3) {
					K[2][j+1] =(pow(s1[2][j+1],2)+pow(s2[2][j+1],2)+pow(s3[2][j+1],2))-(2.0/3.0)*pow(Y[j+1],2);
					if (K[2][j+1] > 0) {
						xx = sqrt(pow(s1[2][j+1],2)+pow(s2[2][j+1],2)+pow(s3[2][j+1],2));
						xx = sqrt(2.0/3.0)*Y[j+1]/xx;
						s1[2][j+1] = xx*s1[2][j+1];
						s2[2][j+1] = xx*s2[2][j+1];
						s3[2][j+1] = xx*s3[2][j+1];
					}
					K[3][j+1] =(pow(s1[3][j+1],2)+pow(s2[3][j+1],2)+pow(s3[3][j+1],2))-(2.0/3.0)*pow(Y[j+1],2);
					if (K[3][j+1] > 0) {
						xx = sqrt(pow(s1[3][j+1],2)+pow(s2[3][j+1],2)+pow(s3[3][j+1],2));
						xx = sqrt(2.0/3.0)*Y[j+1]/xx;
						s1[3][j+1] = xx*s1[3][j+1];
						s2[3][j+1] = xx*s2[3][j+1];
						s3[3][j+1] = xx*s3[3][j+1];
					}
					if(debug == 1) {
						printf("Yield Check %2.0f %2.0f\n", Y[j+1], K[3][j+1]);
					}
				}
				
				//Artificial Viscosity
				//=================================================================================================================
				if((U[2][j+2] < U[2][j]) && (V[3][j+1]-V[1][j+1] < 0)) {
					xx = rho[j+1]/V[2][j+1];
					if (P[1][j+1] < 0) {
						adaptiveMeshRefinement = sqrt(-P[1][j+1]/xx);
					}
					else {
						adaptiveMeshRefinement = sqrt(P[1][j+1]/xx);
					}
					if (fabs(P[1][j+1]) < 0.002) {
						adaptiveMeshRefinement = 0.0;
					}
					q[2][j+1] = artViscConstant1*artViscConstant1*xx*pow((U[2][j+2]-U[2][j]),2)+artViscConstant2*adaptiveMeshRefinement*xx*fabs(U[2][j+2]-U[2][j]);
				}
				else {
					q[2][j+1] = 0.0;
				}
				
				//Energy
				//=================================================================================================================
				
				//ieos(2,imat)=1 - Mie Grunisen
				//ieos(2,imat)=2 - Gamma law ideal gas
				//ieos(2,imat)=3 - Gamma law ideal gas with Newtonian Stress and heat transfer
				//ieos(2,imat)=4 - snow plow model with KO inputs for Hugoniot
				//ieos(2,imat)=5 - snow plow model with anamolous hugoniot
				//ieos(2,imat)=6 - P-alpha Model
				
				if (ieos[1][j+1] == 1) {
					gamma0 = matProp[ieos[0][j+1]][12];
					k1 = rho[j+1]*pow(matProp[ieos[0][j+1]][9],2);
					k2 = k1*(2.0*matProp[ieos[0][j+1]][10]-gamma0/2.0);
					k3 = k1*(3.0*matProp[ieos[0][j+1]][10]-gamma0)*matProp[ieos[0][j+1]][10];
					qbar = (q[2][j+1]+q[0][j+1])/2.0;
					deltaZ = V[2][j+1]*(s1[2][j+1]*epsilon1[2][j+1]+(geometry-1.0)*s2[2][j+1]*epsilon2[2][j+1])*deltat;
					strain = 1.0 - V[3][j+1];
					xb = gamma0;
					if(strain < 0.0) {
						xa = k1*strain+k3*pow(strain,3);
					}
					else {
						xa = k1*strain+k2*pow(strain,2)+k3*pow(strain,3);
					}
					E[3][j+1] = (E[1][j+1]-((xa+P[1][j+1])/2.0+qbar)*(V[3][j+1]-V[1][j+1])+deltaZ)/(1.0+xb*(V[3][j+1]-V[1][j+1])/2.0);
					E[2][j+1] = (E[3][j+1]+E[1][j+1])/2.0;
				}
				else if (ieos[1][j+1] == 2) {
					qbar = (q[2][j+1]+q[0][j+1])/2.0;
					delt = t[2]-t[1];
					deltaZ = V[2][j+1]*(s1[2][j+1]*epsilon1[2][j+1]+(geometry-1.0)*s2[2][j+1]*epsilon2[2][j+1])*delt;
					xa = 0.0;
					xb = (matProp[ieos[0][j+1]][12]-1.0)/V[3][j+1];
					E[3][j+1] = (E[1][j+1]-((xa+P[1][j+1])/2.0+qbar)*(V[3][j+1]-V[1][j+1])+deltaZ)/(1.0+xb*(V[3][j+1]-V[1][j+1])/2.0);
					E[2][j+1] = (E[3][j+1]+E[1][j+1])/2.0;
				}
				else if (ieos[1][j+1] == 3) {
					qbar = (q[2][j+1]+q[0][j+1])/2.0;
					delt = t[2]-t[1];
					deltaZ = V[2][j+1]*(s1[2][j+1]*epsilon1[2][j+1]+(geometry-1.0)*s2[2][j+1]*epsilon2[2][j+1])*delt;
					xx = rho[j+1]/V[3][j+1];
					xa = 0.0;
					xb = (matProp[ieos[0][j+1]][12]-1.0)/V[3][j+1];
					E[3][j+1]= (E[1][j+1]-((xa+P[1][j+1])/2.0+qbar)*(V[3][j+1]-V[1][j+1])+deltaZ)/(1.0+xb*(V[3][j+1]-V[1][j+1])/2.0);
				}
				else if (ieos[1][j+1] == 4) {
					qbar = (q[2][j+1]+q[0][j+1])/2.0;
					deltaZ = V[2][j+1]*(s1[2][j+1]*epsilon1[2][j+1]+(geometry-1.0)*s2[2][j+1]*epsilon2[2][j+1])*deltat;
					xb = matProp[ieos[0][j+1]][12];
					v0 = matProp[ieos[0][j+1]][20];
					
					//Make the compression relative to the compacted material hugoniot
					xx = 1.0-(V[3][j+1]/rho[j+1])/v0;
					if (V[3][j+1]/rho[j+1] <= v0) {
						icompact[3][j+1] = 1;
					}
					if((icompact[3][j+1] == 1) && (xx < 0.0)) {
						xa = matProp[ieos[0][j+1]][9]*xx+matProp[ieos[0][j+1]][11]*pow(xx,3);
					}
					else if ((icompact[3][j+1] == 1) && (xx >= 0.0)) {
						xa = matProp[ieos[0][j+1]][9]*xx+matProp[ieos[0][j+1]][10]*pow(xx,2)+matProp[ieos[0][j+1]][11]*pow(xx,3);
					}
					
					//Zero out pressure until threshold is reached
					else if(icompact[3][j+1] == 0) {
						xa = 0.0;
					}
					else {
						printf("Error - Compaction issue in energy computation\n");
						printf("%d %2.0f\n",icompact[3][j+1], xx);
						errorCode = true;
					}
					E[3][j+1] = (E[1][j+1]-((xa+P[1][j+1])/2.0+qbar)*(V[3][j+1]-V[1][j+1])+deltaZ)/(1.0+xb*(V[3][j+1]-V[1][j+1])/2.0);
					E[2][j+1] = (E[3][j+1]+E[1][j+1])/2.0;
				}
				else if (ieos[1][j+1] == 5) {
					qbar = (q[2][j+1]+q[0][j+1])/2.0;
					deltaZ = V[2][j+1]*(s1[2][j+1]*epsilon1[2][j+1]+(geometry-1.0)*s2[2][j+1]*epsilon2[2][j+1])*deltat;
					stemp = 0.1048;      //slope
					ctemp = 0.5124;      //bulk sound speed
					gtemp = 0.9;         //gamma
					v0 = matProp[ieos[0][j+1]][16];
					v00= 1.0/rho[j+1];
					vv = V[3][j+1]/rho[j+1];    //v_local
					if (vv <= v0) {
						icompact[3][j+1] = 1;
					}
					//icompact is the flag, once compact always compact 
					if (icompact[3][j+1] == 1) { 
						
						//porous hugoniot, see meyers pg 141
						xa = fabs(((2.0*vv-gtemp*(v0-vv))*ctemp*ctemp*(v0-vv))/((2.0*vv-gtemp*(v00-vv))*(pow((v0-stemp*(v0-vv)),2)))); 
					}
					else {
						xa = 0.0;
					}
					E[3][j+1] = (E[1][j+1]-((xa+P[1][j+1])/2.0+qbar)*(V[3][j+1]-V[1][j+1])+deltaZ )/(1.0+xb*(V[3][j+1]-V[1][j+1])/2.0);
					E[2][j+1] = (E[3][j+1] + E[1][j+1])/2.0;
				}
				else if (ieos[1][j+1] == 6) {
					qbar = (q[2][j+1]+q[0][j+1])/2.0;
					deltaZ = V[2][j+1]*(s1[2][j+1]*epsilon1[2][j+1]+(geometry-1.0)*s2[2][j+1]*epsilon2[2][j+1])*deltat;
					xx = 1.0-V[3][j+1];
					xb = matProp[ieos[0][j+1]][12]; //gamma
					if(xx < 0.) {
						xa = matProp[ieos[0][j+1]][9]*xx+matProp[ieos[0][j+1]][11]*pow(xx,3);
					}
					else {
						xa = matProp[ieos[0][j+1]][9]*xx+matProp[ieos[0][j+1]][10]*pow(xx,2)+matProp[ieos[0][j+1]][11]*pow(xx,3);
					}
					E[3][j+1] = (E[1][j+1]-((xa+P[1][j+1])/2.0+qbar)*(V[3][j+1]-V[1][j+1])+deltaZ)/(1.0+xb*(V[3][j+1]-V[1][j+1])/2.0);
					E[2][j+1] = (E[3][j+1]+E[1][j+1])/2.0;
				}
				else {
					printf("Error - EOS computation\n");
					errorCode = true;
				}
				
				//Pressure
				//=================================================================================================================
				
				//Mie Gruneisen EOS as formulated in Wilkins
				if (ieos[1][j+1] == 1) {
					gamma0 = matProp[ieos[0][j+1]][12];
					k1 = rho[j+1]*pow(matProp[ieos[0][j+1]][9],2);
					k2 = k1*(2.0*matProp[ieos[0][j+1]][10]-gamma0/2.0);
					k3 = k1*(3.0*matProp[ieos[0][j+1]][10]-gamma0)*matProp[ieos[0][j+1]][10];
					strain = 1.0-V[2][j+1];

					if (strain < 0.0) {
						P[2][j+1] = k1*strain+k3*pow(strain,3)+gamma0*E[2][j+1];
					}
					else {
						P[2][j+1] = k1*strain+k2*pow(strain,2)+k3*pow(strain,3)+gamma0*E[2][j+1];
					}
					strain = 1.0-V[3][j+1];
					if (strain < 0.0) {
						P[3][j+1]=k1*strain+k3*pow(strain,3)+ gamma0*E[3][j+1];
					}
					else {
						P[3][j+1] = k1*strain+k2*pow(strain,2)+k3*pow(strain,3)+gamma0*E[3][j+1];
					}
				}
				else if (ieos[1][j+1] == 2) {          // Gamma Law (perfect gas)
					P[3][j+1] = (matProp[ieos[0][j+1]+9][12]-1.0)*E[3][j+1]/V[3][j+1];
					P[2][j+1] = (P[1][j+1]+P[3][j+1])/2.0;
				}
				else if (ieos[1][j+1] == 3) {           // Gamma Law
					P[3][j+1] = (matProp[ieos[0][j+1]][12]-1.0)*E[3][j+1]/V[3][j+1];
					P[2][j+1] = (P[1][j+1]+P[3][j+1])/2.0;
				}
				else if (ieos[1][j+1] == 4) {           //Snow Plow
					v0 = matProp[ieos[0][j+1]][20];
					v00 = 1.0/rho[j+1];
					vv = V[2][j+1]/rho[j+1];
					if (v0 == 0) {
						printf("Error - v0 can not equal zero");
						errorCode = true;
					}
					if (vv <= v0) {
						icompact[2][j+1] = 1;
					}
					
					//Make the compression relative to the compacted material hugoniot
					xx = 1.0-(V[2][j+1]/rho[j+1])/v0;
					
					if (icompact[2][j+1] == 1 && xx < 0) {
						P[2][j+1] = matProp[ieos[0][j+1]][9]*xx+matProp[ieos[0][j+1]][11]*pow(xx,3)+matProp[ieos[0][j+1]][12]*E[2][j+1];
					}
					else if(icompact[2][j+1] == 1 && xx >= 0) {
						P[2][j+1] = matProp[ieos[0][j+1]][9]*xx+matProp[ieos[0][j+1]][10]*pow(xx,2)+matProp[ieos[0][j+1]][11]*pow(xx,3)+matProp[ieos[0][j+1]][12]*E[2][j+1];
					}
					else if (icompact[2][j+1] == 0) {
						P[2][j+1] = 0.0*matProp[ieos[0][j+1]][12]*E[2][j+1];
					}
					else {
						printf("Error - Pressure compaction\n");
						errorCode = true;
					}
					if (V[3][j+1]/rho[j+1] <= v0) {
						icompact[3][j+1] = 1;
					}
					
					//Make compression relative to the compacted material hugoniot
					xx = 1.0-(V[3][j+1]/rho[j+1])/v0;
					if (icompact[3][j+1] == 1 && xx < 0) {
						P[3][j+1] = matProp[ieos[0][j+1]][9]*xx+matProp[ieos[0][j+1]][11]*pow(xx,3)+matProp[ieos[0][j+1]][12]*E[3][j+1];
					}
					else if (icompact[3][j+1] == 1 && xx >= 0) {
						P[3][j+1] = matProp[ieos[0][j+1]][9]*xx+matProp[ieos[0][j+1]][10]*pow(xx,2)+matProp[ieos[0][j+1]][11]*pow(xx,3)+matProp[ieos[0][j+1]][12]*E[3][j+1];
					}
					else if (icompact[3][j+1] == 0) { //v .gt. v0
						P[3][j+1] = 0.0*matProp[ieos[0][j+1]][12]*E[3][j+1]; //DH huh ????? *0.0???
					}
					else {
						printf("Error - Pessure compaction\n");
						errorCode = true;
					}
				}
				else if (ieos[1][j+1] == 5) {
					v0 = matProp[ieos[0][j+1]][16];
					v00= 1.0/rho[j+1];
					vv = V[2][j+1]/rho[j+1];    //v_local
					if (vv <= v0) {
						icompact[2][j+1] = 1;
					}
					if (icompact[2][j+1] == 1) {
						
						//porous hugoniot, see meyers pg 141
						xa = ((2.0*vv-gtemp*(v0-vv))*ctemp*ctemp*(v0-vv))/((2.0*vv-gtemp*(v00-vv))*(pow((v0-stemp*(v0-vv)),2)));
						if (xa < 0.0) {
							printf("Compression less than v0 on anamolous hugoniot, 2\n");
							xa = fabs(xa);
						}
					}
					else {
						xa = 0.0;
					}
					
					P[2][j+1] = xa+matProp[ieos[0][j+1]][12]*E[2][j+1];
					
					vv = V[3][j+1]/rho[j+1]; //v_local
					if (vv <= v0) {
						icompact[3][j+1] = 1;
					}
					if (icompact[3][j+1] == 1) {
						//porous hugoniot, see meyers pg 141
						xa = ((2.0*vv-gtemp*(v0-vv))*ctemp*ctemp*(v0-vv))/((2.0*vv-gtemp*(v00-vv))*(pow((v0-stemp*(v0-vv)),2)));
						if (xa < 0.0) {
							printf("Compression less than v0 on anamolous hugoniot, 3\n");
							xa = fabs(xa);
						}
					}
					else {
						xa = 0.0;
					}
				}
				else if (ieos[1][j+1] == 6) { //p-alpha model
					printf("Error - P-Alpha model for pressure has not been implemented\n");
					errorCode = true;
				}
				else if  (ieos[1][j+1] == 7) { //Mie Gruneisen CTH like formulation
					P0 = 0.0;
					E0 = 0.0;
					T0 = 293.0;
					strain = 1.0-V[2][j+1];
					rho_local = rho[j+1]/V[2][j+1];
					up = (U[2][j]+U[2][j+2])/2.0;
					Us = matProp[ieos[0][j+1]][9]+matProp[ieos[0][j+1]][10]*up+matProp[ieos[0][j+1]][11]/matProp[ieos[0][j+1]][9]*pow(up,2);          //local shock speed
					PH = P0+rho_local*Us*up;
					EH = E0+up*up/2.0;
					TH = exp(matProp[ieos[0][j+1]][12]*strain)*T0;
					xx = 1.0-V[3][j+1];
					P[3][j+1] = matProp[ieos[0][j+1]][9]*xx+matProp[ieos[0][j+1]][10]*pow(xx,2)+matProp[ieos[0][j+1]][11]*pow(xx,3)+matProp[ieos[0][j+1]][12]*E[3][j+1];
				}
				else {
					printf("Error - Unexpected condition in pressure computation\n");
					errorCode = true;
				}
			}
			
			//Temperature and Entropy
			//=================================================================================================================
			if (ieos[0][j+1]>1) {
				T[3][j+1] = E[3][j+1]/(rho[j+1]*matProp[ieos[0][j+1]][16]);
				entropy[3][j+1] = entropy[1][j+1]+s1[3][j+1]/T[3][j+1];
			}
			else {
				T[3][j+1] = 0;
				entropy[3][j+1] = 0;
			}
			
			//Pressure and Energy Convergence Check
			//=================================================================================================================
			if (ieos[1][j+1] == 1) {                   //Mie Grunisen
				gamma0 = matProp[ieos[0][j+1]][12];
				k1 = rho[j+1]*pow(matProp[ieos[0][j+1]][9],2);
				k2 = k1*(2.0*matProp[ieos[0][j+1]][10]-gamma0/2.0);
				k3 = k1*(3.0*matProp[ieos[0][j+1]][10]-gamma0)*matProp[ieos[0][j+1]][10];
				qbar = (q[2][j+1]+q[0][j+1])/2.0;
				deltaZ = V[2][j+1]*(s1[2][j+1]*epsilon1[2][j+1]+(geometry-1.0)*s2[2][j+1]*epsilon2[2][j+1])*deltat;
				strain = 1.0-V[3][j+1];
				xb = gamma0; //gamma
				if(strain < 0.0) {
					xa = k1*strain+k3*pow(strain,3);
				}
				else {
					xa = k1*strain+k2*pow(strain,2)+k3*pow(strain,3);
				}
				En2j1 = (E[1][j+1]-((xa+P[1][j+1])/2.0+qbar)*(V[3][j+1]-V[1][j+1])+deltaZ)/(1.0+xb*(V[3][j+1]-V[1][j+1])/2.0);
				diffE = E[3][j+1]-En2j1;
				if (fabs(diffE) > 0) {
					printf("Energy not converged %2.0f %2.0f\n", En2j1, E[3][j+1]);
				}
			}
		}

		//If debug is enabled, output the condition of the debug node
		if (debug == 1) {
			printf("1 %d %2.0f %2.0f %2.0f %2.0f", debugNode, matProp[ieos[0][debugNode+1]][14], V[3][debugNode+1], V[1][debugNode+1], V[2][debugNode+1]);
			printf("s %2.0f %2.0f %2.0f %2.0f", s1[3][debugNode+1], s2[3][debugNode+1], s3[3][debugNode+1], epsilon1[2][debugNode+1]);
			printf("Energy %2.0f %2.0f %2.0f %2.0f", E[3][debugNode+1], E[1][debugNode+1], xx, (2.0*(1.0+xb*(V[3][debugNode+1]-V[1][debugNode+1])/2.0)));
		}
	
		//Spall: Check for tensile fracture (i.e. does the pressure/stress exceed pfrac) and for physical separation
		//=================================================================================================================
		
		//DEV: Shifting the lists with node indexes used here is a loop-dependent operation and cannot be parallelized
		for (int j=0; j<=totalNodes-2; j+=2) {
			if ((-P[3][j+1]+s1[3][j+1] > pfrac[j]) && ibc[j+2] != 2 && ibc[j+2] != 1 && ibc[j-1] != 9) {
				//printf("Before Separating:\n");
				//printf("\tr[3][j-2]=%2.4f;\tU[2][j-2]=%2.4f;\tP[3][j-2]=%2.4f\n", r[3][j-2], U[2][j-2], P[3][j-2]);
				//printf("\tr[3][j-1]=%2.4f;\tU[2][j-1]=%2.4f;\tP[3][j-1]=%2.4f\n", r[3][j-1], U[2][j-1], P[3][j-1]);
				//printf("\tr[3][j]=%2.4f;\tU[2][j]=%2.4f;\tP[3][j]=%2.4f\n", r[3][j], U[2][j], P[3][j]);
				//printf("\tr[3][j+1]=%2.4f;\tU[2][j+1]=%2.4f;\tP[3][j+1]=%2.4f\n", r[3][j+1], U[2][j+1], P[3][j+1]);
				//printf("\tr[3][j+2]=%2.4f;\tU[2][j+2]=%2.4f;\tP[3][j+2]=%2.4f\n", r[3][j+2], U[2][j+2], P[3][j+2]);
				//printf("\tr[3][j+3]=%2.4f;\tU[2][j+3]=%2.4f;\tP[3][j+3]=%2.4f\n", r[3][j+3], U[2][j+3], P[3][j+3]);
				//printf("\tr[3][j+4]=%2.4f;\tU[2][j+4]=%2.4f;\tP[3][j+4]=%2.4f\n", r[3][j+4], U[2][j+4], P[3][j+4]);
				
				for (int jjj=totalNodes-4; jjj>=j+2; jjj-=2) {
					ibc[jjj] = ibc[jjj-2];
					ibc[jjj-1] = ibc[jjj-3];
					ieos[0][jjj] = ieos[0][jjj-2];
					ieos[0][jjj-1] = ieos[0][jjj-3];
					ieos[1][jjj] = ieos[1][jjj-2];
					ieos[1][jjj-1] = ieos[1][jjj-3];
					m[jjj] = m[jjj-2];
					m[jjj-1] = m[jjj-3];
					rho[jjj] = rho[jjj-2];
					rho[jjj-1] = rho[jjj-3];
					Y[jjj] = Y[jjj-2];
					Y[jjj-1] = Y[jjj-3];
					pfrac[jjj] = pfrac[jjj-2];
					pfrac[jjj-1] = pfrac[jjj-3];
				}
				
				#pragma omp parallel for
				for (int nz=1; nz<=3; nz++) {
					for (int jjj=totalNodes-4; jjj>=j+2; jjj-=2) {
						r[nz][jjj] = r[nz][jjj-2];
						r[nz][jjj-1] = r[nz][jjj-3];
						U[nz][jjj] = U[nz][jjj-2];
						U[nz][jjj-1] = U[nz][jjj-3];
						phi[nz][jjj] = phi[nz][jjj-2];
						phi[nz][jjj-1] = phi[nz][jjj-3];
						beta[nz][jjj] = beta[nz][jjj-2];
						beta[nz][jjj-1] = beta[nz][jjj-3];
						sigmar[nz][jjj] = sigmar[nz][jjj-2];
						sigmar[nz][jjj-1] = sigmar[nz][jjj-3];
						sigmao[nz][jjj] = sigmao[nz][jjj-2];
						sigmao[nz][jjj-1] = sigmao[nz][jjj-3];
						V[nz][jjj] = V[nz][jjj-2];
						V[nz][jjj-1] = V[nz][jjj-3];
						s1[nz][jjj] = s1[nz][jjj-2];
						s1[nz][jjj-1] = s1[nz][jjj-3];
						s2[nz][jjj] = s2[nz][jjj-2];
						s2[nz][jjj-1] = s2[nz][jjj-3];
						s3[nz][jjj] = s3[nz][jjj-2];
						s3[nz][jjj-1] = s3[nz][jjj-3];
						q[nz][jjj] = q[nz][jjj-2];
						q[nz][jjj-1] = q[nz][jjj-3];
						epsilon1[nz][jjj] = epsilon1[nz][jjj-2];
						epsilon1[nz][jjj-1] = epsilon1[nz][jjj-3];
						epsilon2[nz][jjj] = epsilon2[nz][jjj-2];
						epsilon2[nz][jjj-1] = epsilon2[nz][jjj-3];
						E[nz][jjj] = E[nz][jjj-2];
						E[nz][jjj-1] = E[nz][jjj-3];
						K[nz][jjj] = K[nz][jjj-2];
						K[nz][jjj-1] = K[nz][jjj-3];
					}
					U[nz][j+1] = 0;
					U[nz][j+1] = 0;
					phi[nz][j+1] = 0;
					beta[nz][j+1] = 0;
					sigmar[nz][j+1] = 0;
					sigmao[nz][j+1] = 0;
					V[nz][j+1] = 0;
					s1[nz][j+1] = 0;
					s2[nz][j+1] = 0;
					s3[nz][j+1] = 0;
					E[nz][j+1] = 0;
					Y[j+1] = 0;
					q[nz][j+1] = 0;
					pfrac[j+1] = 0;
				}
				
				ibc[j] = 2;
				ibc[j+1] = 9;
				ibc[j+2] = -2;
				ieos[0][j+1] = 0;
				ieos[1][j+1] = 0;
				P[2][j+1] = voidPressure;
				P[3][j+1] = voidPressure;
				P[2][j-1] = voidPressure;
				P[3][j-1] = voidPressure;
				P[2][j+3] = voidPressure;
				P[3][j+3] = voidPressure;
				U[2][j+2] = U[2][j+4];
				r[2][j+1] = r[2][j];
				r[3][j+1] = r[3][j];
				pfrac[j] = pfrac[j-1];
				pfrac[j+2] = pfrac[j+3];
				
				//printf("After Separating: %d\n", j);
				//printf("\tr[3][j-2]=%2.4f;\tU[2][j-2]=%2.4f;\tP[3][j-2]=%2.4f\n", r[3][j-2], U[2][j-2], P[3][j-2]);
				//printf("\tr[3][j-1]=%2.4f;\tU[2][j-1]=%2.4f;\tP[3][j-1]=%2.4f\n", r[3][j-1], U[2][j-1], P[3][j-1]);
				//printf("\tr[3][j]=%2.4f;\tU[2][j]=%2.4f;\tP[3][j]=%2.4f\n", r[3][j], U[2][j], P[3][j]);
				//printf("\tr[3][j+1]=%2.4f;\tU[2][j+1]=%2.4f;\tP[3][j+1]=%2.4f\n", r[3][j+1], U[2][j+1], P[3][j+1]);
				//printf("\tr[3][j+2]=%2.4f;\tU[2][j+2]=%2.4f;\tP[3][j+2]=%2.4f\n", r[3][j+2], U[2][j+2], P[3][j+2]);
				//printf("\tr[3][j+3]=%2.4f;\tU[2][j+3]=%2.4f;\tP[3][j+3]=%2.4f\n", r[3][j+3], U[2][j+3], P[3][j+3]);
				//printf("\tr[3][j+4]=%2.4f;\tU[2][j+4]=%2.4f;\tP[3][j+4]=%2.4f\n", r[3][j+4], U[2][j+4], P[3][j+4]);
			}
		}
		
		//Spall End
		dt_min = 1e9;
		#pragma omp parallel for private(Vdot, deltar, b, rho_local, adaptiveMeshRefinement, delt_temp)
		for (int j=0; j<=totalNodes-2; j+=2) {
			if (ibc[j+1] == 0) {
				Vdot = (V[3][j+1]-V[1][j+1])/deltat;
				deltar = fabs(r[3][j+2]-r[3][j]);
				b = 8.0*(pow(artViscConstant1,2)+artViscConstant2)*deltar*(Vdot/V[2][j+1]);
				if (Vdot/V[2][j+1] >= 0) {
					b = 0.0;
				}
				rho_local = rho[j+1]/V[3][j+1];
				if (rho_local > 1e10) {
					printf("Error: Density too high\n");
					printf("%d %2.4f\n",j,t[1]);
					printf("%2.4f %2.4f\n", rho[j+1], V[3][j+1]);
					printf("%2.4f %2.4f\n", r[3][j+2], r[3][j]);
					printf("%2.4f %2.4f\n", V[3][j+1], V[1][j+1]);
					printf("%2.4f %2.4f %2.4f %2.4f\n", P[2][j], P[2][j+1], P[3][j], P[3][j+1]);
					errorCode = true;
				}
				adaptiveMeshRefinement = sqrt(fabs(P[3][j+1])/rho_local);
				if ((pow(adaptiveMeshRefinement,2)+pow(b,2)) != 0) {
					delt_temp = (2.0/3.0)*(deltar/sqrt(pow(adaptiveMeshRefinement,2)+pow(b,2)));
					#pragma omp critical
					{
					if (delt_temp < dt_min) {
						dt_min = delt_temp;
					}
					}
				}
			}
		}
		
		
		//Adjust true timestep according to the adaptive mesh refinement value
		if (deltat > dt_min/6) {
			deltat = dt_min/6;
		}
		if (deltat < dt_min/20) {
			deltat = dt_min/20;
		}
		
		//Output solution to console
		//=================================================================================================================
		if ((int)((-1)/printIterSkip) == -1/printIterSkip || t[3] >= tstop) {
			qtotal  = 0.0;  // total artifical viscosity
			mvtotal = 0.0;  // total momentum
			etotal  = 0.0;  // total energy
			ketotal = 0.0;  // total kinetic energy
			ietotal = 0.0;  // total internal energy
			ke3total= 0.0;  // total kinetic energy
			
			#pragma omp parallel for reduction(+:qtotal) reduction(+:mvtotal) reduction(+:ketotal) reduction(+:etotal) reduction(+:ke3total)
			for (int j=0; j<=totalNodes-2; j+=2) {
				if (ibc[j+1] == 0) {
					qtotal += q[3][j];
					mvtotal += m[j+1]*(U[2][j]+U[2][j+2])/2.0;
					ketotal += m[j+1]*pow(((U[2][j]+U[2][j+2])/2.0),2)/2.0;
					etotal += (E[3][j+1]-E[1][j+1])/deltat+P[2][j+1]*(V[3][j+1]-V[1][j+1])/deltat;
					if (ieos[0][j] == 2) {
						ke3total += m[j+1]*pow(((U[2][j]+U[2][j+2])/2.0),2)/2.0;
					}
				}
			}
			printf("%d %2.15f %2.15f %2.15f %2.15f %2.15f %2.15f\n", timestep, t[2], mvtotal, ketotal, ietotal, etotal, deltat);
		}
		
		
		//Write solution to file
		//=================================================================================================================
		if (t[1] == 0 || t[1] >= fileWrite) {
			//DEV: Not sure of the exact file output format, need to update accordingly before parallelizing
			for (int j=0; j<=totalNodes-4; j+=2) {
				bs1 = V[1][j+1]; //rho(j+1)
				bs2 = U[2][j];
				bs3 = P[1][j+1];
				bs4 = E[1][j+1];
				bs5 = q[2][j+1];
				bs6 = s1[1][j+1];
				bs7 = s2[1][j+1];
				bs8 = epsilon1[0][j+1];
				bs10 = T[1][j+1];
				
				//If the exponent is larger than 99, zero it out
				if (fabs(bs1) < 1e-99) {
					bs1 = 0.0;
				} 
				if (fabs(U[2][j]) < 1e-99) {
					bs2 = 0.0;
				}
				if (fabs(bs3) < 1e-99) {
					bs3 = 0.0;
				}
				if (fabs(bs4) < 1e-99) {
					bs4 = 0.0;
				}
				if (fabs(bs5) < 1e-99) {
					bs5 = 0.0;
				}
				if (fabs(bs6) < 1e-99) {
					bs6 = 0.0;
				}
				if (fabs(bs7) < 1e-99) {
					bs7 = 0.0;
				}
				if (fabs(bs8) < 1e-99) {
					bs8 = 0.0;
				}
				bs9 = -bs3+bs6;
			}
			fileWrite = fileWrite+fileTimeSkip;
		}
		
		//Update solution
		//=================================================================================================================
		#pragma omp parallel for
		for (int j=0; j<totalNodes; j++) {
			if (ibc[j] != 9) {
				U[0][j] = U[2][j];
				U[1][j] = U[3][j];
				phi[0][j] = phi[2][j];
				phi[1][j] = phi[3][j];
				sigmar[0][j] = sigmar[2][j];
				sigmar[1][j] = sigmar[3][j];
				sigmao[0][j] = sigmao[2][j];
				sigmao[1][j] = sigmao[3][j];
				beta[0][j] = beta[2][j];
				beta[1][j] = beta[3][j];
				V[0][j] = V[2][j];
				V[1][j] = V[3][j];
				r[0][j] = r[2][j];
				r[1][j] = r[3][j];
				epsilon1[0][j] = epsilon1[2][j];
				epsilon1[1][j] = epsilon1[3][j];
				epsilon2[0][j] = epsilon2[2][j];
				epsilon2[1][j] = epsilon2[3][j];
				s1[0][j] = s1[2][j];
				s1[1][j] = s1[3][j];
				s2[0][j] = s2[2][j];
				s2[1][j] = s2[3][j];
				s3[0][j] = s3[2][j];
				s3[1][j] = s3[3][j];
				P[0][j] = P[2][j];
				P[1][j] = P[3][j];
				q[0][j] = q[2][j];
				q[1][j] = q[3][j];
				t[0] = t[2];
				t[1] = t[3];
				E[0][j] = E[2][j];
				E[1][j] = E[3][j];
				K[0][j] = K[2][j];
				K[1][j] = K[3][j];
				T[0][j] = T[2][j];
				T[1][j] = T[3][j];
				entropy[0][j] = entropy[2][j];
				entropy[1][j] = entropy[3][j];
			}
		}
		
		//Check for the occurance of an error generated during execution of the parallelized sections
		if (errorCode) {
			printf("Error - An error occured during a parallelized section during the last loop\n");
			return -1;
		}
		
	}
	
	//Output the complete program computation time
	gettimeofday(&endTime, 0);
	printf("Computation time (s): %f\n", (((endTime.tv_sec-startTime.tv_sec)*1000.0)+((endTime.tv_usec-startTime.tv_usec)/1000.0))/1000);
	
	//Terminate program
	return 0;
	
}
//=================================================================================================================
