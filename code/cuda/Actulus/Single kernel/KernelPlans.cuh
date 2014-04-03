#include <device_launch_parameters.h>
__global__ void RunPlans(float *peVa, float *peResult, 
						 float *dtlaVa, float *dtlaResult,
						 float *tlapVa, float *tlapResult, 
						 float *tiVa, float *tiResult, 
						 float *daVa, float *daResult, 
						 float *dtiVa, float *dtiResult);