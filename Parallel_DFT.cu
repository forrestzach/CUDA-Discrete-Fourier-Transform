    //Forrest Zach
	//CPSC 4780
	//Final Project
	
	// the code calculates a DFT of a random complex number input and 
    // then an IDFT. The IDFT result should be the input vector 
    // to compile with gcc
    // gcc -Wall -O2 -fopenmp -o DFTW DFTW.c 
    // written by stef

    // exercise


	#include <stdio.h> // printf
	#include <stdlib.h> // malloc and rand for instance. Rand not thread safe!
	#include <sys/time.h>   // time(0) to get random seed
	#include "math.h"  // sine and cosine
	// #include "omp.h"   // openmp library like timing
	#include <cuda_runtime.h>

	
	// two pi
	#define PI2 6.28318530718
	// this for the rounding error, increasing N rounding error increases
    // 0.01 precision good for N > 8000
	#define R_ERROR 0.01
	
	// set the input array with random number
	int fillInput(double* xr, double* xi, int N, int whatInput);
	// set to zero the input vector
	int setOutputZero(double* Xr_o, double* Xi_o, int N);
	// check if x = IDFT(DFT(x))
	int checkResults(double* xr, double* xi, double* xr_check, double* xi_check, double* Xr_o, double* Xi_r, int N, int isCPUCall);
	// print the results of the DFT
	int printResults(double* xr, double* xi, int N, int isCPUCall);

	// main routine to calculate DFT
	// DFT/IDFT routine
	// idft: 1 direct DFT, -1 inverse IDFT (Inverse DFT)
	int DFT(int idft, double* xr, double* xi, double* Xr_o, double* Xi_o, int N){
	  for (int k=0 ; k<N ; k++)
	    {
	        for (int n=0 ; n<N ; n++)  {
	        	// Real part of X[k]
	            Xr_o[k] += xr[n] * cos(n * k * PI2 / N) + idft*xi[n]*sin(n * k * PI2 / N);
	            // Imaginary part of X[k]
	            Xi_o[k] += -idft*xr[n] * sin(n * k * PI2 / N) + xi[n] * cos(n * k * PI2 / N);
	            
	        } 
	    }
	    
	    // normalize if you are doing IDFT
	    if (idft==-1){
	    	for (int n=0 ; n<N ; n++){
	    	Xr_o[n] /=N;
	    	Xi_o[n] /=N;
	    }
	    }
	  return 1; 
	}
	// 16 * 512 = 8192, I believe that is the amount of threads that will be computing here
	__global__ void dft_kernel(int idft, double* xr, double* xi, double* Xr_o, double* Xi_o, int N){
		unsigned int tid = threadIdx.x;
    	unsigned int idx = blockIdx.x * blockDim.x + tid;

		// printf("blockIdx.x: %d, blockDim.x: %d, threadIdx.x: %d, idx: %d \n", blockIdx.x, blockDim.x, tid, idx);
		//so I have idx which can be the index of each part of the output
		//now with that, I need to iterate through the input data and run the calculations
		//basically just the nested for loop only, referring to output[idx]

	    if(idx >= N) return; //The relevance of this line of code makes much more sense

		for (int n=0 ; n<N ; n++)  {
	       	// Real part of X[idx]
			Xr_o[idx] += xr[n] * cos(n * idx * PI2 / N) + idft*xi[n]*sin(n * idx * PI2 / N);
			// Imaginary part of X[idx]
			Xi_o[idx] += -idft*xr[n] * sin(n * idx * PI2 / N) + xi[n] * cos(n * idx * PI2 / N);
			
		}
		//Normalize part for IDFT 
		if (idft==-1){
			// printf("In normalize part\n");
	    	Xr_o[idx] /=N;
			// printf("Xr_o[idx]: %lf, Xr_o[idx] / N: %lf\n", Xr_o[idx], Xr_o[idx]/N);
	    	Xi_o[idx] /=N;
	    }
	}


	int main(int argc, char* argv[]){
	  // size of input array
	  int N = 8000; // 8,000 is a good number for testing
	  if(!(argc < 2)){
		  N = atoi(argv[1]);
	  }
	  else{
		  printf("Must enter N value! Command should look like ./a.out <N Value> <OPTIONAL Input Choice>\n");
		  return 0;
	  }
	  
	  printf("DFTW calculation with N = %d \n",N);
	  int whatInput = 1;
	  if(argc > 2){
		  whatInput = atoi(argv[2]);
	  }

	  // execution configuration
      int blocksize = 512;   // initial block size

      dim3 block (blocksize, 1);
      dim3 grid  ((N + block.x - 1) / block.x, 1);
      printf("grid %d block %d\n", grid.x, block.x);
	  
	  // Allocate array for input vector
	  double* xr = (double*) malloc (N *sizeof(double));
	  double* xi = (double*) malloc (N *sizeof(double));
	  fillInput(xr,xi,N, whatInput);

	  // for checking purposes
	  double* xr_check = (double*) malloc (N *sizeof(double));
	  double* xi_check = (double*) malloc (N *sizeof(double));
	  setOutputZero(xr_check,xi_check,N);
	  
	  // Allocate array for output vector
	  double* Xr_o = (double*) malloc (N *sizeof(double));
	  double* Xi_o = (double*) malloc (N *sizeof(double));
	  setOutputZero(Xr_o,Xi_o,N);
	  
	  //Prepare device data
	  size_t bytes = N * sizeof(double);

      cudaEvent_t start, stop;
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
	  //Establish and allocate device memory
	  double *d_xr = NULL;
	  double *d_xi = NULL;
	  double *d_xr_check = NULL;
	  double *d_xi_check = NULL;
	  double *d_Xr_o = NULL;
	  double *d_Xi_o = NULL;

	  cudaMalloc((void **) &d_xr, bytes);
	  cudaMalloc((void **) &d_xi, bytes);
	  cudaMalloc((void **) &d_xr_check, bytes);
	  cudaMalloc((void **) &d_xi_check, bytes);
	  cudaMalloc((void **) &d_Xr_o, bytes);
	  cudaMalloc((void **) &d_Xi_o, bytes);

	  //I recognize it is strange to do all variable copying here, but
	  //	I want the vectors to be zeroed out
	  
	  cudaMemcpy(d_xr, xr, bytes, cudaMemcpyHostToDevice);
	  cudaMemcpy(d_xi, xi, bytes, cudaMemcpyHostToDevice);
	  cudaMemcpy(d_xr_check, xr_check, bytes, cudaMemcpyHostToDevice);
	  cudaMemcpy(d_xi_check, xi_check, bytes, cudaMemcpyHostToDevice);
	  cudaMemcpy(d_Xr_o, Xr_o, bytes, cudaMemcpyHostToDevice);
	  cudaMemcpy(d_Xi_o, Xi_o, bytes, cudaMemcpyHostToDevice);

	  //CPU Computations

	  int idft = 1;
	  //CPU timer
	  struct timeval timestamp;
  	  gettimeofday(&timestamp, NULL);
  	  double timer1 = (double)timestamp.tv_sec + ((double)timestamp.tv_usec / 100000.0);
	  // DFT
	  DFT(idft,xr,xi,Xr_o,Xi_o,N);
	  // IDFT
	  idft = -1;
	  DFT(idft,Xr_o,Xi_o,xr_check,xi_check,N);

      gettimeofday(&timestamp, NULL);
      double timer2 = (double)timestamp.tv_sec + ((double)timestamp.tv_usec / 100000.0);
      double differenceTime = timer2-timer1;

	  // check the results: easy to make correctness errors with openMP
	  checkResults(xr,xi,xr_check,xi_check,Xr_o, Xi_o, N, 1);
	  // print the results of the DFT
      printResults(Xr_o,Xi_o,N, 1);

	  //Clear output data from CPU run
	  setOutputZero(xr_check,xi_check,N);
	  setOutputZero(Xr_o,Xi_o,N);


	  //GPU Computations
	  cudaDeviceSynchronize();
	  idft = 1;

	  cudaEventRecord(start);
	  //DFT Kernel Call
	  dft_kernel<<<grid, block>>>(idft, d_xr, d_xi, d_Xr_o, d_Xi_o, N);

	  cudaDeviceSynchronize();
	  idft = -1;
	  //IDFT Kernel Call
	  dft_kernel<<<grid, block>>>(idft, d_Xr_o, d_Xi_o, d_xr_check, d_xi_check, N);
	
	  cudaEventRecord(stop);

	  cudaDeviceSynchronize();
	  //Copy back GPU Output
	  cudaMemcpy(Xr_o, d_Xr_o, bytes, cudaMemcpyDeviceToHost);
	  cudaMemcpy(Xi_o, d_Xi_o, bytes, cudaMemcpyDeviceToHost);
	  cudaMemcpy(xr_check, d_xr_check, bytes, cudaMemcpyDeviceToHost);
	  cudaMemcpy(xi_check, d_xi_check, bytes, cudaMemcpyDeviceToHost);
	  //Perform check on GPU output to ensure that it is correct
	  checkResults(xr,xi,xr_check,xi_check,Xr_o, Xi_o, N, 0);

	//   printf("GPU Xre[1000] = %f \n",Xr_o[1000]);

      float milliseconds = 0.0;
      cudaEventElapsedTime(&milliseconds, start, stop);
	  printResults(Xr_o,Xi_o,N,0);

	  printf("\nTotal CPU processing time(s): %f\n", differenceTime);
      printf("Total GPU processing time(s): %f\n", (milliseconds/1000));

	  printf("GPU Speedup over CPU: %f\n", differenceTime/(milliseconds/1000));
	  
	  // take out the garbage
	  free(xr); free(xi);
	  free(Xi_o); free(Xr_o);
	  free(xr_check); free(xi_check);

	  return 1;
	}
	
	// set the initial signal 
    // be careful with this 
    // rand() is NOT thread safe in case
	//whatInput determines which input generation to use, 1 (default) uses random data
	//	2 uses the constant real signal
	int fillInput(double* xr, double* xi, int N, int whatInput){
	    srand(time(0));
	    for(int n=0; n < 100000;n++) // get some random number first 
	    	rand();
	    for(int n=0; n < N;n++){
	       // Generate random discrete-time signal x in range (-1,+1)
		   if(whatInput == 1){
	         xr[n] = ((double)(2.0 * rand()) / RAND_MAX) - 1.0;
	         xi[n] = ((double)(2.0 * rand()) / RAND_MAX) - 1.0;
		   }
	       // constant real signal
		   if(whatInput == 2){
	       	 xr[n] = 1.0;
	         xi[n] = 0.0;
		   }
	    }
		return 1; 
	}

	// set to zero the output vector
	int setOutputZero(double* Xr_o, double* Xi_o, int N){
	for(int n=0; n < N;n++){
	       Xr_o[n] = 0.0;
	       Xi_o[n] = 0.0; 
	    }
		return 1;
	}

	// check if x = IDFT(DFT(x))
	int checkResults(double* xr, double* xi, double* xr_check, double* xi_check, double* Xr_o, double* Xi_r, int N, int isCPUCall){
		// x[0] and x[1] have typical rounding error problem
		// interesting there might be a theorem on this
		for(int n=0; n < N;n++){
			if (fabs(xr[n] - xr_check[n]) > R_ERROR)
			    printf("ERROR - x[%d] = %f, inv(X)[%d]=%f \n",n,xr[n], n,xr_check[n]);
			if (fabs(xi[n] - xi_check[n]) > R_ERROR)
			    printf("ERROR - x[%d] = %f, inv(X)[%d]=%f \n",n,xi[n], n,xi_check[n]);

		}
		if(isCPUCall){
		    // printf("\nCPU Xre[0] = %f \n",Xr_o[0]);
		}
		else{
			// printf("\nGPU Xre[0] = %f \n",Xr_o[0]);
		}
		return 1;
	}

	// print the results of the DFT
	int printResults(double* xr, double* xi, int N, int isCPUCall){
		if(isCPUCall){
			printf("\nBegin CPU Data!\n");
			printf("Xre[%d] = %f, Xim[%d] = %f \n", 0, xr[0], 0, xi[0]);
		}
		else{
			printf("\nBegin GPU Data!\n");
			printf("Xre[%d] = %f, Xim[%d] = %f \n", 0, xr[0], 0, xi[0]);
		}
		for(int n=N-20; n < N;n++)
			    printf("Xre[%d] = %f, Xim[%d] = %f \n", n, xr[n], n, xi[n]);
		return 1;
	}