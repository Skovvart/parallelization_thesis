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
					printf("  %20.12f ", result[offset + (i * states) + s]);
				}
				printf("\n");
			}
			printf("\n");
		}
	}
}

void runKernel(bool printresult){
	int blocks = 1;
	int threads = 1;
	// Allocate result memory on host
	float *peResult   =	(float*) malloc(sizeof(float) * 41     * blocks * threads);
	float *dtlaResult = (float*) malloc(sizeof(float) * 51     * blocks * threads);
	float *tlapResult = (float*) malloc(sizeof(float) * 51     * blocks * threads);
	float *tiResult   =	(float*) malloc(sizeof(float) * 51	   * blocks * threads);
	float *daResult   =	(float*) malloc(sizeof(float) * 51 * 2 * blocks * threads);
	float *dtiResult  = (float*) malloc(sizeof(float) * 51 * 2 * blocks * threads);

	//Allocate Va on host
	float *peVa		= (float *) malloc(sizeof(float));
	float *dtlaVa	= (float *) malloc(sizeof(float));
	float *tlapVa	= (float *) malloc(sizeof(float));
	float *tiVa		= (float *) malloc(sizeof(float));
	float *daVa		= (float *) malloc(sizeof(float) * 2);
	float *dtiVa	= (float *) malloc(sizeof(float) * 2);
	peVa[0] = 0.0f;
	dtlaVa[0] = 0.0f;
	tlapVa[0] = 0.0f;
	tiVa[0] = 0.0f;
	daVa[0] = 0.0f;
	daVa[1] = 0.0f;
	dtiVa[0] = 0.0f;
	dtiVa[1] = 0.0f;

	// Allocate memory on device, and copy over initial data
	float *d_peResult, *d_dtlaResult, *d_tlapResult, *d_tiResult, *d_daResult, *d_dtiResult;
	float *d_peVa, *d_dtlaVa, *d_tlapVa, *d_tiVa, *d_daVa, *d_dtiVa;

	HANDLE_ERROR(cudaMalloc((void**)&d_peResult,	sizeof(float) * 41 * blocks * threads));
	HANDLE_ERROR(cudaMalloc((void**)&d_dtlaResult,	sizeof(float) * 51 * blocks * threads));
	HANDLE_ERROR(cudaMalloc((void**)&d_tlapResult,	sizeof(float) * 51 * blocks * threads));
	HANDLE_ERROR(cudaMalloc((void**)&d_tiResult,	sizeof(float) * 51 * blocks * threads));
	HANDLE_ERROR(cudaMalloc((void**)&d_daResult,	sizeof(float) * 51 * 2 * blocks * threads));
	HANDLE_ERROR(cudaMalloc((void**)&d_dtiResult,	sizeof(float) * 51 * 2 * blocks * threads));

	HANDLE_ERROR(cudaMalloc((void**)&d_peVa,	sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&d_dtlaVa,	sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&d_tlapVa,	sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&d_tiVa,	sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&d_daVa,	sizeof(float) * 2));
	HANDLE_ERROR(cudaMalloc((void**)&d_dtiVa,	sizeof(float) * 2));

	HANDLE_ERROR(cudaMemcpy(d_peVa,		peVa,	sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_dtlaVa,	dtlaVa,	sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_tlapVa,	tlapVa,	sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_tiVa,		tiVa,	sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_daVa,		daVa,	sizeof(float) * 2, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_dtiVa,	dtiVa,	sizeof(float) * 2, cudaMemcpyHostToDevice));

	// Start the kernel
	timerstart();
	RunPlans<<<blocks, threads>>>(
		d_peVa, d_peResult, 
		d_dtlaVa, d_dtlaResult, 
		d_tlapVa, d_tlapResult, 
		d_tiVa, d_tiResult, 
		d_daVa, d_daResult, 
		d_dtiVa, d_dtiResult);
	kernelTime += timerstop();
	HANDLE_ERROR(cudaPeekAtLastError());
	HANDLE_ERROR(cudaThreadSynchronize());
	//HANDLE_ERROR(cudaGetLastError());


	// Copy the results to main memory
	HANDLE_ERROR(cudaMemcpy(peResult, d_peResult,		sizeof(float) * 41 * blocks * threads, cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(dtlaResult, d_dtlaResult,	sizeof(float) * 51 * blocks * threads, cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(tlapResult, d_tlapResult,	sizeof(float) * 51 * blocks * threads, cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(tiResult, d_tiResult,		sizeof(float) * 51 * blocks * threads, cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(daResult, d_daResult,		sizeof(float) * 51 * 2 * blocks * threads, cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(dtiResult, d_dtiResult,		sizeof(float) * 51 * 2 * blocks * threads, cudaMemcpyDeviceToHost));
	//HANDLE_ERROR(cudaMemcpy(Va, d_Va, sizeof(float) * states, cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaFree(d_peResult));
	HANDLE_ERROR(cudaFree(d_dtlaResult));
	HANDLE_ERROR(cudaFree(d_tlapResult));
	HANDLE_ERROR(cudaFree(d_tiResult));
	HANDLE_ERROR(cudaFree(d_daResult));
	HANDLE_ERROR(cudaFree(d_dtiResult));

	HANDLE_ERROR(cudaFree(d_peVa));
	HANDLE_ERROR(cudaFree(d_dtlaVa));
	HANDLE_ERROR(cudaFree(d_tlapVa));
	HANDLE_ERROR(cudaFree(d_tiVa));
	HANDLE_ERROR(cudaFree(d_daVa));
	HANDLE_ERROR(cudaFree(d_dtiVa));

	if(printresult){
		print(blocks, threads, 40, 1, peResult);
		print(blocks, threads, 50, 1, dtlaResult);
		print(blocks, threads, 50, 1, tlapResult);
		print(blocks, threads, 50, 1, tiResult);
		print(blocks, threads, 50, 2, daResult);
		print(blocks, threads, 50, 2, dtiResult);
	}

	free(peResult);
	free(dtlaResult);
	free(tlapResult);
	free(tiResult);
	free(daResult);
	free(dtiResult);

	free(peVa);
	free(dtlaVa);
	free(tlapVa);
	free(tiVa);
	free(daVa);
	free(dtiVa);
}

int main(int argc, const char* argv[]){
	cudainit();
	runKernel(true);
	printf("Total kernel-time %fms\n", kernelTime);
	return 0;
}