#include <stdio.h>
#include <cuda_runtime.h>

cudaEvent_t start, stop;
float timertime;

void timerstart() {
	timertime = 0.0f;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
}

float timerstop() {
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&timertime, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	return timertime;
}