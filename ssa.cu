#include"1dIsing.h"


/*

	SSA algorithm running on one processor. 

*/ 



double randU(void);
float computeRatesSSA(int*lat, int index, int N,  reactData react);




float ssa(int* lat, int N, int finalT, reactData react){

int i, j, k;
float* rates = NULL;
float u, t;
float rateSum = 0.0;
float tempSum;
int cov = 0;

	rates = (float*)malloc(N * sizeof(float));
	if(rates == NULL){
		printf("ERROR : ssa : No memory for rates\n");
		exit(1);
	}	
		
	for(i = 0; i < N;i++){
		rates[i] = computeRatesSSA(lat, i, N, react);
		rateSum += rates[i];
	}	
	
	t = -log(randU())/rateSum;
	
	for(i = 0;i < N;i++)
			cov += lat[i];
	
	//for(j = 0;j < finalT;j++){
	while(t < finalT){
	
		tempSum = 0.0;
		u = randU() * rateSum;
		
		for(i = 0;i < N; i++){
		
			tempSum += rates[i];
			if(tempSum > u){
			
				lat[i]  = 1 - lat[i];
				rates[i] = computeRatesSSA(lat, i, N, react);
				
				cov += 2*lat[i] - 1;
				
				if(i == (N-1)){
					rates[N - 2] = computeRatesSSA(lat, N - 2 , N, react);
					rates[0]     = computeRatesSSA(lat, 0, N, react);
				}
				else if(i == 0){
					rates[1]     = computeRatesSSA(lat, 1, N, react);
					rates[N - 1] = computeRatesSSA(lat, N - 1, N, react); 
				}
				else{
					rates[i - 1] = computeRatesSSA(lat, i - 1, N, react);
					rates[i + 1] = computeRatesSSA(lat, i + 1, N, react);
				}
				
				
				// Recalculate sum of rates
				rateSum = 0.0;
				for(k = 0; k < N; k++){
					rateSum += rates[k];
				}
				
				break;
			}
		
	   }
	   	   	
		t += -log(randU())/rateSum;
	
	   
	}
	
	
	
	   
	   
	
	free(rates);	
	return cov/(float)(N);
 
}



//===================================================================

float computeRatesSSA(int*Lat, int index, int N,  reactData react){
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

//====================================================================

double randU(void){

	return (double)rand()/(RAND_MAX + 1.0);

}



