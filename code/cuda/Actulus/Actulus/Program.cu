#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include "CudaTimer.cuh"
#include "KernelPlans.cuh"
float kernelTime = 0.0;


static void HandleError(cudaError_t res, const char *file, int line) {
	if (res != cudaSuccess) {
		printf("%s in %s at line %d\n", cudaGetErrorString(res), file, line);
		exit(0);
	}
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

void cudainit() {
	cudaDeviceProp prop;
	int count, usedDevice;
	HANDLE_ERROR(cudaGetDeviceCount(&count));

	int lastclock = 0;
	for (int i=0; i< count; i++) {
		HANDLE_ERROR(cudaGetDeviceProperties(&prop, i));
		// Choose the fastest card
		if(prop.clockRate > lastclock) {
			HANDLE_ERROR(cudaSetDevice(i));
			lastclock = prop.clockRate;
		}    
	}

	HANDLE_ERROR(cudaGetDevice(&usedDevice));
	HANDLE_ERROR(cudaGetDeviceProperties(&prop, usedDevice));
	//printf("Device %d: %s - %dKHz, %dmp\n", usedDevice, prop.name, prop.clockRate, prop.multiProcessorCount);
}

void print(int blocks, int threads, int years, int states, float *result) {	
	for(int block = 0; block < blocks; block++) {
		for(int thread = 0; thread < threads; thread++) {
			unsigned long offset = (thread + block * threads) * states * (years+1);
			for(int i = 0; i <= years; i++) {
				printf("%3d:", i);
				for(int s = 0; s < states; s++) {
					printf("  %20.16f", result[offset + (i * states) + s]);
				}
				printf("\n");
			}
			printf("\n");
		}
	}
}

void runKernel(void (*plan)(int a, int b, int steps, float *Va, float *result), int states, int years, int steps, int blocks, int threads, bool printresult){
	// Calculate the size needed to store results
	unsigned long int resultSize = sizeof(float) * (years+1) * states * blocks * threads;
	// Allocate result memory on host
	float *result = (float*) malloc(resultSize);

	//Allocate Va on host
	float *Va = (float *) malloc(sizeof(float) * states);
	int i;
	for (i = 0; i < states; i++)
		Va[i] = 0.0f;
	
	// Allocate memory on device, and copy over initial data
	float *d_result;
	float *d_Va;

	HANDLE_ERROR(cudaMalloc((void**)&d_result, resultSize));
	HANDLE_ERROR(cudaMalloc((void**)&d_Va, sizeof(float) * states));
	HANDLE_ERROR(cudaMemcpy(d_Va, Va, sizeof(float) * states, cudaMemcpyHostToDevice));

	// Start the kernel
	timerstart();
	plan<<<blocks, threads>>>(years, 0, steps, d_Va, d_result);
	kernelTime += timerstop();
	HANDLE_ERROR(cudaPeekAtLastError());

	// Copy the results to main memory
	HANDLE_ERROR(cudaMemcpy(result, d_result, resultSize, cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(Va, d_Va, sizeof(float) * states, cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaFree(d_result));
	HANDLE_ERROR(cudaFree(d_Va));
	
	if(printresult)
		print(blocks, threads, years, states, result);

	free(result);
	free(Va);
}

void runKernels(bool printresult){
	int steps = 100;
	int blocks = 1;
	int threads = 1;

	
	runKernel(&PureEndowment,					1, 40, steps, blocks, threads, printresult);
	runKernel(&DeferredTemporaryLifeAnnuity,	1, 50, steps, blocks, threads, printresult);
	runKernel(&TemporaryLifeAnnuityPremium,		1, 50, steps, blocks, threads, printresult);
	runKernel(&TermInsurance,					1, 50, steps, blocks, threads, printresult);
	runKernel(&DisabilityAnnuity,				2, 50, steps, blocks, threads, printresult);
	runKernel(&DisabilityTermInsurance,			2, 50, steps, blocks, threads, printresult);

}

int main(int argc, const char* argv[]){
	cudainit();
	int iterations = 10;
	for (int i = 0; i < iterations; i++)
		runKernels(i==0);
	printf("Total kernel-time %fms\n", kernelTime/iterations);
	return 0;
}