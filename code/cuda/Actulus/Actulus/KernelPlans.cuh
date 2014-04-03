#include <device_launch_parameters.h>
__global__ void PureEndowment(int a, int b, int steps, float* Va, float* result);
__global__ void DeferredTemporaryLifeAnnuity(int a, int b, int steps, float* Va, float* result);
__global__ void TemporaryLifeAnnuityPremium(int a, int b, int steps, float* Va, float* result);
__global__ void TermInsurance(int a, int b, int steps, float* Va, float* result);
__global__ void DisabilityAnnuity(int a, int b, int steps, float* Va, float* result);
__global__ void DisabilityTermInsurance(int a, int b, int steps, float* Va, float* result);