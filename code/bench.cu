#include"1dIsing.h"


/* Benchmarks */


void benchH(int* lat, int* devLat, int N,float deltaT, 
		int *cov, int*devCov, reactData react, FILE* results, curandState *devStates, arguments args){

	int j = 0;
	
	int steps = (int)(args.finalT/deltaT);
	float meanCov;
	float ssaCov;
	int *ssaLat = (int*)malloc(N * sizeof(int));
	
	
	
	if(ssaLat == NULL){
		printf("ssaLat = NULL\n");
		exit(1);
	}
		
	for(j = 0;j < N;j++)
		ssaLat[j] = lat[j];
	
	printf("\n    h\t\t DevCov\t\t SSACov\t\tExact solution\tError\n");
	
	j = 0;	
	for(j = 0; j < STEPS_IN_H; j++){

		react.h = j*2.0/(STEPS_IN_H - 1);
	
		//Burning time
		 fractionalStep(lat,devLat,BURNING, react, devStates, deltaT, N, cov, devCov, args);
		 
		//Actual computation
		 meanCov = fractionalStep(lat,devLat,steps,
							react, devStates, deltaT, N, cov, devCov,args);
		
		// ssaCov = ssa(ssaLat, N, 10000, react); 

		fprintf(results,"%f\t%f\n",react.h,meanCov);
		printf("%f\t%f\t%f\t%f\t%1.2e\n",react.h,meanCov,ssaCov,
			exactMagnetization(react),fabs(meanCov-exactMagnetization(react)));
	}

	free(ssaLat);
}

//============================================================


void benchError(int* lat, int* devLat, int N,float deltaT, 
		int *cov, int*devCov, reactData react, FILE* results, curandState *devStates, arguments args){

	int i, j;
	int top = 10;  // Change the top variable to take more points in the benchmark
	float ssaCov;
	int nReals = 20; // Number of realizations to take
	
	int steps = (int)(args.finalT/deltaT);
	
	float meanCov;
	
	if(N!= (int)pow(2,top)){
		printf("ERROR : benchError :\n\tDimension of lattice : %d\n\tDimension required for this test : %d\n\t",N, (int)pow(2,top));
		printf("Either change the 'top' variable in 'benchError'\n\tor the size of the lattice\n\t");
		printf("In the second case, split %d between threads (-t) and blocks (-b)\n",(int)pow(2,top)/8);
		exit(1);
	}
	
	if(args.blocks != 2){
		printf("ERROR : benchError :\tThe number of blocks must be two for this benchmark.\n");
		exit(1);
	}



	ssaCov = ssa(lat, N, 10000, react); 
	printf("Offset\t <m> (fs)\t <m> (ssa)\t exact\t\t Error\n");
	for(i = 2;i < (top - 1);i++){
		args.offSet   = (int)pow(2, i );
		args.threads = (int)pow(2, top - i - 2);
		args.blocks  = 2;
		meanCov = 0.0;
		
		for(j = 0;j < nReals;j++){
			//Burning time
			 //fractionalStep(lat,devLat,BURNING, react, devStates, deltaT, N, cov, devCov, args);
		 
			//Actual computation
			 meanCov += fractionalStep(lat,devLat,steps, react, devStates, deltaT, N, cov, devCov,args);
		}

	    meanCov = meanCov/nReals;
		
		fprintf(results,"%d\t\t%f\n",args.offSet, meanCov, ssaCov,fabs(meanCov - ssaCov));
		printf("%d\t%f\t%f\t%f\t%1.2e\n",args.offSet,meanCov,ssaCov, exactMagnetization(react),fabs(meanCov-exactMagnetization(react)));
		 
	}


}

//===========================================================
/* 

	Benchmark for different DT and constant h

*/


void benchDt(int* lat, int* devLat, int N,float deltaT,
                int *cov, int*devCov, reactData react, FILE* results, curandState *devStates, arguments args){

	int j,i;	
	float meanCov;
	float ssaCov;
	float exactSolution;
	int *ssaLat = (int*)malloc(N * sizeof(int));
	

	const int steps = (int)(args.finalT/deltaT);
	const int DTsteps      = 20;  // Different Dt to try 
	const float largestDt  = 3.0; 
	const float smallestDt = 1/2.0; 
	const int nReals  = 10; // Number of realizations to take
	

	if(ssaLat == NULL){
		printf("ssaLat = NULL\n");
		exit(1);
	}
		
	for(j = 0;j < N;j++)
		ssaLat[j] = lat[j];
	
	printf("\n    Dt\t\t DevCov\t\t SSACov\t\tExact solution\tError\n");
	
	
	react.h = 2*2.0/(STEPS_IN_H - 1);	

	ssaCov = ssa(ssaLat, N, 10000, react); 
	exactSolution = exactMagnetization(react);

	for(j = DTsteps; j > 0; j--){

		args.deltaT = j * (largestDt - smallestDt)/DTsteps + smallestDt;
		meanCov = 0.0;
		for(i = 0; i < nReals; i++){
			//Burning time
			 fractionalStep(lat,devLat,BURNING, react, devStates, deltaT, N, cov, devCov, args);
		 
			//Actual computation
		 	meanCov += fractionalStep(lat,devLat,steps,
							react, devStates, deltaT, N, cov, devCov,args);
		}
		 
		meanCov = meanCov/nReals;
		
		fprintf(results,"%f\t%f\n",args.deltaT,meanCov,exactSolution,fabs(meanCov-exactSolution));
		printf("%f\t%f\t%f\t%f\t%1.2e\n",args.deltaT, meanCov, ssaCov,
			exactSolution,fabs(meanCov-exactSolution));
	}

	free(ssaLat);




}

//============================================================

/* 
	Comparison between exact solution and approximate solution while changing both
    Dt and h (potential).

    */


void benchDt2(int* lat, int* devLat, int N,float deltaT,
                int *cov, int*devCov, reactData react, FILE* results, curandState *devStates, arguments args){


	int j,i;	
	float meanCov;
	float exactSolution;
	

	const int steps = (int)(args.finalT/deltaT);
	const int DTsteps      = 3;  // Different Dt to try 
	const float largestDt  = 3.0; 
	const float smallestDt = 1/2.0; 
	
	
	printf("\n    Dt\t\tDevCov\t\tExact solution\tError\n");
	
	
	for(j = DTsteps; j > 0; j--){

		args.deltaT = j * (largestDt - smallestDt)/DTsteps + smallestDt;
		
		for(i = 0; i < STEPS_IN_H; i++){
			react.h = i * 2.0/(STEPS_IN_H-1);
			//Burning time
			 fractionalStep(lat,devLat,BURNING, react, devStates, deltaT, N, cov, devCov, args);
		 
			//Actual computation
		 	meanCov = fractionalStep(lat,devLat,steps,
							react, devStates, deltaT, N, cov, devCov,args);

		 	exactSolution = exactMagnetization(react);
		 	fprintf(results,"%f\t%f\t%f\t%f\t%f\n",args.deltaT,react.h,meanCov,exactSolution,fabs(meanCov-exactSolution));
		    printf("%f\t%f\t%f\t%1.2e\n",args.deltaT, meanCov,exactSolution,fabs(meanCov-exactSolution));
		}
		 
	}

}



//============================================================

//% exact solution of the mean coverage of the Ising model
//% with spins in {0,1}
float exactMagnetization(reactData react){
float b = react.b;
float h = react.h;
float K01 = react.J;
h = -h/2+K01/2;
float K = K01/4;

float a = sinh(b*h);
      b = a*a + exp(-b*4*K);
float y = a/sqrt(b);
y = (y+1)/2;
return y;
}



