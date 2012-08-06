#include"1dIsing.h"

/*
	Main device code for the 1D fractional step
*/



//=========================================================================

//Amazing helper function. Use it to get the correct errors back from the 
//CUDA functions.
#define checkCudaErrors(err)  __checkCudaErrors (err, __FILE__, __LINE__)

inline void __checkCudaErrors(cudaError err, const char *file, const int line )
{
    if(cudaSuccess != err)
    {
        fprintf(stderr, "%s(%i) : CUDA Runtime API error %d: %s.\n",file,
							line, (int)err, cudaGetErrorString( err ) );
        exit(-1);        
    }
}
//========================================================================







int main(int argc, char* argv[]){

arguments args;

parser(argc, argv, &args); // parses command line arguments

int N = 2 * args.offSet * args.blocks * args.threads; // N : Size of the system
int NTHREADS = args.threads;
int NBLOCKS  = args.blocks;


int *lat = (int*)malloc(N * sizeof(int)); // local lattice
int *cov; //= (int*)malloc(NBLOCKS*sizeof(int));

int *devLat, i, *devCov;
float time;

cudaEvent_t start, stop;
//States for the random number generator
curandState *devStates;



reactData react;


//cov = (int*)malloc(NBLOCKS*sizeof(int));

// Allocate pinned memory for better performance 


cudaMallocHost(&cov, args.blocks * sizeof(int));	
checkCudaErrors(cudaMalloc((void**)&devCov, args.blocks * sizeof(int)));
zeroCov<<<args.blocks, args.threads>>>(devCov); // Initialize devCov variable with zeroes


//==============================================

FILE *results = NULL;
	
	
	results = fopen("RESULTS.txt","w");
	if(results == NULL){
		printf("Couldn't open file for saving the results. Will exit\n");	
		exit(-1);
	}

//===============================================
	//Initialize the react structure
	react.h  = 0.0;
	react.b  = 2.0;
	react.J  = 1.0;
	react.cd = 1.0;
	react.ca = 1.0; 
//==============================================

if(argc < 2)
	printf(" This program can take arguments from the command line.\n Use --help to see the available options.\n");


float deltaT = args.deltaT;

//=====================================================

//Check if there is an actual device in the system

int cudaCount;
checkCudaErrors(cudaGetDeviceCount(&cudaCount));

if(cudaCount == 0){
	printf(" Program can't run without a GPU. Will now exit.\n");
	exit(1);
}
else{
	printf(" Found at least one device!\n\n");
}

//====================================================


// Allocate space for the lattice on the device
checkCudaErrors(cudaMalloc((void**)&devLat, N * sizeof(int)));


// Allocate space for the RNG states on the device
checkCudaErrors(cudaMalloc((void**)&devStates, args.threads * args.blocks * sizeof(curandState)));


//=====================================================

printf("============================================================\n\t\t\tFractional step 1D\n============================================================\n\n\n");
printf(" Final time we want to reach : %d\n", args.finalT);
printf(" Deltat : %1.2f\n",args.deltaT);
printf(" Steps to reach Final Time : %d\n", (int)(args.finalT/args.deltaT));
printf(" Blocks used  : %d\n",NBLOCKS);
printf(" Threads used : %d\n", NTHREADS);
printf(" Number of cells in each part : %d\n", args.offSet);
printf(" System size    : %d\n", N);
printf(" Splitting used : %s\n", args.splitName);

//========================================

//Initialize the lattice with values
for(i = 0; i < N; i++){
	lat[i] = 1;
} 


//=======================================


printf(" Lattice initialized. Locked & ready for computation!\n");

printf("==============================================\n\n");

// Copy the lattice onto the GPU
		checkCudaErrors(cudaMemcpy(devLat, lat, N * sizeof(int), cudaMemcpyHostToDevice));

// Initalize the random number generators for the threads	
		initRandomSequences<<<NBLOCKS,NTHREADS>>>(devStates,1403);
		
		
//Initiate events and measure execution time
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);	

// For benchmarks look in : bench.cu

	//benchH(lat, devLat, N, deltaT, cov,  devCov, react, results, devStates, args);

	//benchError(lat, devLat, N, deltaT, cov,  devCov, react, results, devStates, args);
	//benchDt(lat, devLat, N, deltaT, cov,  devCov, react, results, devStates, args);
	benchDt2(lat, devLat, N, deltaT, cov,  devCov, react, results, devStates, args);
    
    cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	time = time/(1000*60); // Translate from millisec to minutes
	printf("\n Execution time : %f mins\n",time);



//======================================
//Clean on the card side
printf(" Freeing memory on the card side . . . ");
cudaFree(devLat);
cudaFree(devCov);
printf("Success\n");
//=======================================

//Clean on the host side 
printf(" Freeing memory on the host side . . . ");
free(lat);

// Free pinned memory
cudaFreeHost(cov);
printf("Success\n");
//=======================================

printf(" Closing files . . . ");
//Close the file
fclose(results);
printf("Success\n");

printf("================ Program finished ==============\n\n");

return 0;
}



