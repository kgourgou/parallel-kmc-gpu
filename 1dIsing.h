#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<curand.h>
#include<curand_kernel.h>


#define DtoH cudaMemcpyDeviceToHost

// How many points to be used for the external field. Used for the mean coverage benchmark. 
#define STEPS_IN_H 20

#define BURNING 30


//=========================================================


// enum variable for the different splittings. You need to add new splittings here.
enum splittings{

	Lie,
	Strang,
	Random

};

//Struct that holds all the data for a particular reaction
typedef struct
{
	float h;
	float b;
	float J;
	float cd;
	float ca;

} reactData;

// Struct that holds data from the command line
typedef struct {

int threads;
int blocks;
int offSet;
int finalT;
float deltaT;
enum splittings split;
char splitName[10];

} arguments;




//============ Host functions ====================================================

float fractionalStep(int*Lat, int* devLat,int steps, reactData react,
	 		curandState *devStates, float deltaT, int N, int* cov, int*devCov, arguments args);
float exactMagnetization(reactData react);
void parser(int argc, char* argv[], arguments* args);

//================ Benchmarks : for explanations look in bench.cu

void benchH(int* lat, int* devLat, int N,float deltaT,int *cov, int*devCov, 
		reactData react, FILE*  results,curandState *devStates, arguments args);
void benchError(int* lat, int* devLat, int N,float deltaT, 
		int *cov, int*devCov, reactData react, FILE* results, curandState *devStates,arguments args);
void benchDt(int* lat, int* devLat, int N,float deltaT,
                int *cov, int*devCov, reactData react, FILE* results, curandState *devStates, arguments args);
void benchDt2(int* lat, int* devLat, int N,float deltaT,
                int *cov, int*devCov, reactData react, FILE* results, curandState *devStates, arguments args);

//================ Splittings
void chooseSplit(int* devLat, reactData react, curandState *devStates, int*devCov, arguments args);
void splitLie   (int* devLat, reactData react, curandState *devStates, int*devCov, arguments args);
void splitStrang(int* devLat, reactData react, curandState *devStates, int*devCov, arguments args);
void splitRandom(int* devLat, reactData react, curandState *devStates, int*devCov, arguments args);

//================ SSA ===========================================================
float ssa(int* lat, int N, int finalT, reactData react);

//=========== Device functions ===================================================
__global__ void L(int* Lat, int N, float Dt, int Number_Of_Partitions, int var, reactData react,curandState *globalState, int*devCov);
__device__ float computeRates(int* Lat, int index, int N, reactData react);
__global__ void initRandomSequences(curandState *state, int seed);
__global__ void calcCov(int *Lat, int offSet ,int *devCov);
__global__ void zeroCov(int* devCov);

