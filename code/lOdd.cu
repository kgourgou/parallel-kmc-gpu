#include"1dIsing.h"


//N here is the system size constant


//var = 0 : odd
//var = 1 : even

__global__ void L(int* Lat, int N, float Dt, int NumOfPar, int var,reactData react ,curandState *globalState, int* devCov){

int localIndex = NumOfPar * (2 * threadIdx.x + var) + NumOfPar * 2 * blockDim.x * blockIdx.x;

int idofX = threadIdx.x + blockIdx.x * blockDim.x;
int totalCov = 0;
int i;
float t = 0.0, rateSum, u, temp;


curandState s  = globalState[idofX];

//Sum of rates
		for(i = localIndex; i < (localIndex + NumOfPar);i++){
				rateSum  += computeRates(Lat, i, N, react);
				totalCov += Lat[i];
			}
	
		t = -log(curand_uniform(&s))/rateSum;
	
	while(t < Dt){
	
		//Do KMC on that part
		temp = 0.0;
		u =  curand_uniform(&s);
		u = rateSum * u;
		
		for(i = localIndex; i < (localIndex + NumOfPar); i++){

			temp  += computeRates(Lat,i,N,react);
			
			if(temp > u){
			
				Lat[i]    = 1 - Lat[i];	       // spin - flip	
				totalCov += (2 * Lat[i]) - 1; // update the local coverage
				break;
			}
		}
		
		// Sample from the exponential distribution and add to the total time
		rateSum = 0.0;
		for(i = localIndex; i < (localIndex + NumOfPar);i++){
				rateSum += computeRates(Lat, i, N, react);
			}
			
		t += -log(curand_uniform(&s))/rateSum;
	}


	

        globalState[idofX] = s;
}

//=========================================================================



/*

 Reduces the computation of the coverage by first computing the each (thread, sublattice)
 and then each (block, sublattice). The CPU then just has to sum a 1 x #blocks array instead
 of a 1 x N lattice. 
 
 */
__global__ void calcCov(int *Lat, int offSet ,int *devCov){
int localIndexEven = offSet * (2 * threadIdx.x + 0) + offSet * 2 * blockDim.x * blockIdx.x;
int localIndexOdd  = offSet * (2 * threadIdx.x + 1) + offSet * 2 * blockDim.x * blockIdx.x;
int j;
extern __shared__ int buffer[]; 
int totalCov = 0;
	
	
	for(j = localIndexEven; j < (localIndexEven + offSet); j++){ // Even indexes of the lattice for each thread
				totalCov += Lat[j];
			}
			
	for(j = localIndexOdd; j < (localIndexOdd + offSet);j++){ // Odd indexes of the lattice for each thread
		totalCov += Lat[j];	
			}
	       
		buffer[threadIdx.x] = totalCov;  // Save each local coverage to shared memory
	 	__syncthreads();
	 	
	 	totalCov = 0;
		if(threadIdx.x == 0){ // Calculate the coverage per block
		
			for(j = 0; j < blockDim.x; j++){		
		        	totalCov += buffer[j];
			}
		
			devCov[blockIdx.x]  += totalCov; // Save the coverage in global memory
		}
}

//=========================================================================
__global__ void zeroCov(int* devCov){

	if(threadIdx.x == 0)
		devCov[blockIdx.x] = 0;
}



//=========================================================================
// Computes the rates. This version works only for nearest neighbor
// interactions.



__device__ float computeRates(int* Lat, int index, int N, reactData react){
int prIndex, neIndex;

//==================================
// Chemical reaction specific constants
float h  = react.h;
float b  = react.b;
float J  = react.J;
float cd = react.cd;
float ca = react.ca;
//==================================

float U,r;

	if(index == 0 ){
		prIndex = (N - 1);
		neIndex = 1;
	}
	else if(index == (N - 1)){
		prIndex = (N - 2);
		neIndex = 0;
	}
	else{
		prIndex = (index - 1);
		neIndex = (index + 1);
	}


U = J * (Lat[prIndex] + Lat[neIndex]) - h;

r = cd *exp(-b * U) * Lat[index] + ca * (1 - Lat[index]) ;

return  r;
}


//================================================================

// Setup the random number generators
__global__ void initRandomSequences(curandState *state, int seed)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;

		curand_init(seed, id , 0, &state[id]) ;

	}



