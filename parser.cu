
#include"1dIsing.h"
#include<getopt.h>
#include <string.h>

/*

	parser : parses the arguments from the command line and saves the results in the 
	"arguments" structure.

*/
void parser(int argc, char* argv[], arguments* args){

int optionIndex = 0;
int c  = 0 ;
char* shortOpts = "t:,b:,o:,T:,D:,h:,s:"; // Short options 
int choice;

//========================================================

struct option longOpts[] ={ // Long options
 {"threads", required_argument, 0, 't'},
 {"blocks" , required_argument, 0, 'b'},
 {"offSet" , required_argument, 0, 'o'},
 {"time"   , required_argument, 0, 'T'},
 {"deltaT" , required_argument, 0, 'd'},
 {"help"   , no_argument      , 0, 'h'},
 {"split"  ,required_argument , 0, 's'},
 {   0	   ,         0        , 0,  0}};

//======================= Help  text =========================

 char* help = "============================================================\n\t\t\tFractional step 1D\n============================================================\n\nAvailable options : \n\n--threads = Number of threads per block to use (short : -t)\n--blocks  = Number of blocks to use (short -b)\n--offSet  = Size of sublattice to use per thread (short -o)\n--time    = Final time to reach (short -T)\n--deltaT  = deltaT with which to progress (short --d)\n\nExample:\n\n ./a.out -t 10 --d 2\n will run the program with 10 threads,\n 1 block, 4x4 sublattices and a Dt of 2.0."; 
 

//==================== Default values ==============================

	args->threads = 10;
	args->blocks  = 1;
	args->offSet  = 4;
	args->finalT  = 1000;
	args->deltaT  = 1.0;
	args->split   = Lie; 


//===================================================


while(1){
	c = getopt_long(argc, argv, shortOpts, longOpts, &optionIndex);
	if(c == -1)  break;
	switch (c)
		     {
		
		     case 't': // Number of threads
		       args->threads = atoi(optarg);
		       break;
	     
		     case 'b': // Number of blocks
			   args->blocks = atoi(optarg);
		       break;
	     
		     case 'o': // Offset used
		       args->offSet = atoi(optarg);
		       break;
	     
		     case 'T': // Final Time
		       args->finalT = atoi(optarg);
		       break;
	     
		     case 'd': // DeltaT 
		       args->deltaT = atof(optarg);
		       break;

		     case 's': // Splitting. You need to add new splittings here

		        choice = atoi(optarg);
		        if(choice == Lie){
		        	args->split = Lie;
		        	strcpy(args->splitName, "Lie");
		        }
		        else if(choice == Strang){
		        	args->split = Strang;
		        }
		        else if(choice == Random){
		        	args->split = Random;
		        }
		        else{
		        	printf("ERROR : parser : unrecognized value given for the splitting\n");
		        	printf("Program will now exit\n");
		        	exit(1);
		        }
		     	break;

		     case 'h':// help 
		     	puts(help);
		     	exit(0);
	         	break;
	             
		     default:
				printf("Found unrecognized option\n");
				exit(1);
				break;
		     }
		 }
}

