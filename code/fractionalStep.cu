#include"1dIsing.h"


/* 

	Fractional step algorithm. Look near the end of the file for the different splittings available.

*/

float fractionalStep(int* lat, int* devLat,int steps, reactData react,
			 curandState *devStates, float deltaT, int N, int *cov, int*devCov, arguments args){
int i,l;
float meanCov = 0.0;
int newCoverage = 0; 
size_t sharedMem = args.threads * sizeof(int); // Shared memory size



for(i = 0; i < steps; i++){

			chooseSplit(devLat,react,devStates,devCov,args);
	       	
			// Calculate the mean coverage of the last (STEPS_IN_TIME - BURNING) samples	
	        calcCov<<<args.blocks, args.threads, sharedMem>>>(devLat, args.offSet, devCov);                
	}	
	        	        
	/*cudaMemcpy(lat, devLat, N*sizeof(int),DtoH);
	for(l = 0;l<N;l++){
		newCoverage += lat[l];
		}
		
		
		printf("newCoverage = %d\n", newCoverage);		
	*/	
	cudaMemcpy(cov, devCov, args.blocks * sizeof(int),DtoH);
	
	for(l = 0; l < args.blocks; l++){
		newCoverage += cov[l];
	}	


       meanCov = newCoverage/(N*(steps)*1.0); 
       zeroCov<<<args.blocks, args.threads>>>(devCov);
	
       
return meanCov;
	
}

//========================================================================
/*

	Splittings : Lie, Strang, Random	
	
	All splittings have the same argument list, so only the change of the name is needed.
	
*/


void chooseSplit(int* devLat, reactData react, curandState *devStates, int *devCov,
												 arguments args){

	if(args.split == Lie){
 		splitLie(devLat, react, devStates, devCov, args);
	}
	else if(args.split == Strang)
	{
		splitStrang(devLat, react, devStates, devCov, args);
	}
	else if(args.split == Random)
	{

		splitRandom(devLat, react, devStates, devCov, args);

	}
}

//================================

void splitLie(int* devLat, reactData react, curandState *devStates, int *devCov, arguments args){
	int nthreads = args.threads;
	int nblocks  = args.blocks;
	int offSet   = args.offSet;
	int N = 2 * args.offSet * args.blocks * args.threads; // N : Size of the system
	float deltaT = args.deltaT;
	
      			L<<<nblocks, nthreads>>>(devLat, N, deltaT, offSet,0, react,devStates,devCov);
		        L<<<nblocks, nthreads>>>(devLat, N, deltaT, offSet,1, react,devStates,devCov);
		
}
//=================================

void splitStrang(int* devLat, reactData react, curandState *devStates, int *devCov, arguments args){

	int nthreads = args.threads;
	int nblocks  = args.blocks;
	int offSet   = args.offSet;
	int N = 2 * args.offSet * args.blocks * args.threads; // N : Size of the system
	float deltaT = args.deltaT;


			L<<<nblocks, nthreads>>>(devLat, N, deltaT/2.0, offSet, 1,react, devStates, devCov);
		  L<<<nblocks, nthreads>>>(devLat, N, deltaT,     offSet, 0,react, devStates, devCov);
			L<<<nblocks, nthreads>>>(devLat, N, deltaT/2.0, offSet, 1,react, devStates, devCov);

}

//================================
void splitRandom(int* devLat, reactData react, curandState *devStates, int *devCov, arguments args){

	int nthreads = args.threads;
	int nblocks  = args.blocks;
	int offSet   = args.offSet;
	int N = 2 * args.offSet * args.blocks * args.threads; // N : Size of the system
	float deltaT = args.deltaT;

	
	if(rand() < (1/2.0)){
		
		L<<<nblocks, nthreads>>>(devLat, N, deltaT, offSet,0,react,devStates, devCov);
		
	}
	else{
	
		L<<<nblocks, nthreads>>>(devLat, N, deltaT, offSet,1, react,devStates, devCov);
		
	}
}


//===============================













