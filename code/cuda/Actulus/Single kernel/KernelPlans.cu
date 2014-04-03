#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

__device__ float age = 30.0f;
__device__ float interestrate = 0.05f;
__device__ float bpension = 1.0f;
__device__ float pensiontime = 35.0f;

__device__ float indicator(int b) {
	return b ? 1.0f : 0.0f;
}

// Gompertz-Makeham mortality intensities for Danish women
__device__ float GM(float t) {
	return 0.0005f + pow(10.0f, 5.728f - 10.0f + 0.038f * (age + t));
}

__device__ float r(float t) {
	return interestrate;    // Fixed interest rate
	// return rFsa(t);            // Finanstilsynet's rate curve
}

// The Danish FSA yield curve (Finanstilsynets rentekurve).
// Data from 2011-11-16 
__device__ float ts[] = { 
	0.25f,  0.5f,  1.0f,  2.0f,  3.0f,  4.0f,  5.0f,  6.0f,  7.0f,  8.0f, 9.0f, 10.0f, 
	11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 16.0f, 17.0f, 18.0f, 19.0f, 20.0f,
	21.0f, 22.0f, 23.0f, 24.0f, 25.0f, 26.0f, 27.0f, 28.0f, 29.0f, 30.0f 
};

__device__ float rs[] = { 
	1.146677033f, 1.146677033f, 1.146677033f, 1.340669678f, 1.571952911f, 1.803236144f, 
	2.034519377f, 2.265802610f, 2.497085843f, 2.584085843f, 2.710085843f, 2.805085843f, 
	2.871485843f, 2.937885843f, 3.004285843f, 3.070685843f, 3.137085843f, 3.136485843f, 
	3.135885843f, 3.135285843f, 3.134685843f, 3.134085843f, 3.113185843f, 3.092285843f, 
	3.071385843f, 3.050485843f, 3.029585843f, 3.008685843f, 2.987785843f, 2.966885843f, 
	2.945985843f, 2.925085843f
};

__device__ float rFsa(float t){
	// Requires ts non-empty and elements strictly increasing.
	int last = sizeof(ts)/sizeof(float) - 1;
	if(t <= ts[0])
		return log(1 + rs[0] / 100);
	if(t >= ts[last])
		return log(1 + rs[last] / 100);
	int a = 0; 
	int b = last;
	// Now a < b (bcs. ts must have more than 1 element) and ts[a] < t < ts[b]
	while(a + 1 < b) {
		// Now a < b and ts[a] <= t < ts[b]
		int i = (a + b) / 2;
		if(ts[i] <= t)
			a = i;
		else // t < ts[i]
			b = i;
	}
	// Now a+1>=b and ts[a] <= t < ts[b]; so a!=b and hence a+1 == b <= last
	int m = a;
	float tm = ts[m];
	float tm1 = ts[m + 1];
	float rm = rs[m] / 100; 
	float rm1 = rs[m + 1] / 100;
	float Rt = (rm * (tm1 - t) + rm1 * (t - tm)) / (tm1 - tm);
	return log(1 + Rt) + t / (tm1 - tm) * (rm1 - rm) / (1 + Rt);
}

// sax = scalar a times x array, imperative version
__device__ void sax(float a, float* x, float* res, int states) {
	for(int i = 0; i < states; i++, x++)
		res[i] = a * *x;
}

// saxpy = scalar a times x array plus y array, imperative version
__device__ void saxpy(float a, float* x, float* y, float* res, int states) {
	for(int i = 0; i < states; i++, x++, y++)
		res[i] = a * *x + *y;
}

__device__ void floatArrayCopy(float* from, float* to, int length){
	for(int i = 0; i < length; i++)
		to[i] = from[i];
}

__device__ unsigned long getOffset(int thread, int block, int blocks, int states, int year){
	return (thread + block * blocks) * states * (year + 1);
}

//begin PureEndowment
#define PureEndowmentStates 1

__device__ float pe_b_0(float t) { return 0.0f; }
__device__ float pe_mu_01(float t) { return GM(t); }
__device__ float pe_bj_00(float t) { return t == pensiontime ? bpension : 0.0f; }
__device__ float pe_bj_01(float t) { return 0.0f; }

__device__ void pe_dV(float t, float* V, float* res){
	res[0] = r(t) * V[0] - pe_b_0(t) - pe_mu_01(t) * (0 - V[0] + pe_bj_01(t));
}

__device__ void pe_bj_ii(float t, float* res){
	res[0] = pe_bj_00(t);
}
//end PureEndowment

//begin DeferredTemporaryLifeAnnuity
#define DeferredTemporaryLifeAnnuityStates 1
__device__ int dtla_m = 35;
__device__ int dtla_n = 10;
__device__ float dtla_b_0(float t) { return bpension * indicator(t > dtla_m) * indicator(t < dtla_m + dtla_n); }
__device__ float dtla_mu_01(float t) { return GM(t); }
__device__ float dtla_bj_00(float t) { return 0.0f; }
__device__ float dtla_bj_01(float t) { return 0.0f; }
__device__ void dtla_dV(float t, float* V, float* res){ res[0] = r(t) * V[0] - dtla_b_0(t) - dtla_mu_01(t) * (0 - V[0] + dtla_bj_01(t)); }
__device__ void dtla_bj_ii(float t, float* res){ res[0] = dtla_bj_00(t); }
//end DeferredTemporaryLifeAnnuity

//begin TemporaryLifeAnnuityPremium
#define TemporaryLifeAnnuityPremiumStates 1
__device__ int tlap_n = 35;
__device__ float tlap_bpremium = 1.0f;

__device__ float tlap_b_0(float t) { return -tlap_bpremium * indicator(t >= 0) * indicator(t < tlap_n); }
__device__ float tlap_mu_01(float t) { return GM(t); }
__device__ float tlap_bj_00(float t) { return 0.0f; }
__device__ float tlap_bj_01(float t) { return 0.0f; }

__device__ void tlap_dV(float t, float* V, float* res){
	res[0] = r(t) * V[0] - tlap_b_0(t) - tlap_mu_01(t) * (0 - V[0] + tlap_bj_01(t));
}

__device__ void tlap_bj_ii(float t, float* res){
	res[0] = tlap_bj_00(t);
}
//end TemporaryLifeAnnuityPremium

//begin TermInsurance
#define TermInsuranceStates 1
__device__ int ti_n = 35;
__device__ float ti_bdeath = 1.0f;

__device__ float ti_b_0(float t) { return 0.0f; }
__device__ float ti_mu_01(float t) { return GM(t); }
__device__ float ti_bj_00(float t) { return 0.0f; }
__device__ float ti_bj_01(float t) { return ti_bdeath * indicator(t > 0) * indicator(t < ti_n); }

__device__ void ti_dV(float t, float* V, float* res){
	res[0] = r(t) * V[0] - ti_b_0(t) - ti_mu_01(t) * (0 - V[0] + ti_bj_01(t));
}

__device__ void ti_bj_ii(float t, float* res){
	res[0] = ti_bj_00(t);
}
//end TermInsurance

//begin DisabilityAnnuity
#define DisabilityAnnuityStates 2
__device__ int da_n = 35;
__device__ float da_bdisabled = 1.0f;

__device__ float da_b_0(float t) { return 0.0f; }
__device__ float da_b_1(float t) { return da_bdisabled * indicator(t > 0) * indicator(t < da_n); }
__device__ float da_GM01(float t) { return 0.0006f + pow(10.0f, 4.71609f - 10.0f + 0.06f * (age + t)); }
__device__ float da_GM02(float t) { return GM(t); }
__device__ float da_GM12(float t) { return GM(t); }
__device__ float da_mu_01(float t) { return da_GM01(t); }
__device__ float da_mu_02(float t) { return da_GM02(t); }
__device__ float da_mu_12(float t) { return da_GM12(t); }
__device__ float da_bj_00(float t) { return 0.0f; }
__device__ float da_bj_01(float t) { return 0.0f; }
__device__ float da_bj_02(float t) { return 0.0f; }
__device__ float da_bj_11(float t) { return 0.0f; }
__device__ float da_bj_12(float t) { return 0.0f; }

__device__ void da_dV(float t, float* V, float* res){
	res[0] = r(t) * V[0] - da_b_0(t) - da_mu_01(t) * (V[1] - V[0] + da_bj_01(t)) - da_mu_02(t) * (0 - V[0] + da_bj_02(t));
	res[1] = r(t) * V[1] - da_b_1(t) - da_mu_12(t) * (0 - V[1] + da_bj_12(t)); 
}

__device__ void da_bj_ii(float t, float* res){
	res[0] = da_bj_00(t); 
	res[1] = da_bj_11(t);
}
//end DisabilityAnnuity

//begin DisabilityTermInsurance
#define DisabilityTermInsuranceStates 2
__device__ int dti_n = 35;
__device__ float dti_bdisabled = 1.0f;

__device__ float dti_b_0(float t) { return 0.0f; }
__device__ float dti_b_1(float t) { return 0.0f; }
__device__ float dti_GM01(float t) { return 0.0006 + pow(10.0f, 4.71609f - 10.0f + 0.06f * (age + t)); }
__device__ float dti_GM02(float t) { return GM(t); }
__device__ float dti_GM12(float t) { return GM(t); }
__device__ float dti_mu_01(float t) { return dti_GM01(t); }
__device__ float dti_mu_02(float t) { return dti_GM02(t); }
__device__ float dti_mu_12(float t) { return dti_GM12(t); }
__device__ float dti_bj_00(float t) { return 0.0f; }
__device__ float dti_bj_01(float t) { return dti_bdisabled * indicator(t > 0) * indicator(t < dti_n); }
__device__ float dti_bj_02(float t) { return 0.0f; }
__device__ float dti_bj_11(float t) { return 0.0f; }
__device__ float dti_bj_12(float t) { return 0.0f; }

__device__ void dti_dV(float t, float* V, float* res){
	res[0] = r(t) * V[0] - dti_b_0(t) - dti_mu_01(t) * (V[1] - V[0] + dti_bj_01(t)) - dti_mu_02(t) * (0 - V[0] + dti_bj_02(t));	
	res[1] = r(t) * V[1] - dti_b_1(t) - dti_mu_12(t) * (0 - V[1] + dti_bj_12(t));
}

__device__ void dti_bj_ii(float t, float* res){
	res[0] = dti_bj_00(t); 
	res[1] = dti_bj_11(t);
}
//end DisabilityTermInsurance

__device__ void PureEndowment(int a, int b, int steps, float* Va, float* result){
	unsigned long offset = getOffset(threadIdx.x, blockIdx.x, blockDim.x, PureEndowmentStates, a);

	float h = -1.0f / steps;
	float k1[PureEndowmentStates], k2[PureEndowmentStates], k3[PureEndowmentStates], k4[PureEndowmentStates], tmp[PureEndowmentStates], v[PureEndowmentStates];

	floatArrayCopy(Va, &result[offset + a * PureEndowmentStates], PureEndowmentStates);

	for(int y = a; y > b; y--) {
		pe_bj_ii(y, v);
		saxpy(1.0, v, &result[offset + y * PureEndowmentStates], v, PureEndowmentStates);
		for(int s = 0; s < steps; s++) { 	// Integrate backwards over [y,y-1]
			float t = y - s / (float)steps;
			// Hack: Fake limit from left
			pe_dV(s == 0 ? t - 1e-5 : t, v, k1);
			sax(h, k1, k1, PureEndowmentStates);
			saxpy(0.5, k1, v, tmp, PureEndowmentStates);
			pe_dV(t + h / 2, tmp, k2);
			sax(h, k2, k2, PureEndowmentStates);
			saxpy(0.5, k2, v, tmp, PureEndowmentStates);
			pe_dV(t + h / 2, tmp, k3);
			sax(h, k3, k3, PureEndowmentStates);
			saxpy(1, k3, v, tmp, PureEndowmentStates);
			// Hack: Fake limit from right
			pe_dV(s == steps - 1 ? t + h + 1e-5 : t + h, tmp, k4);
			sax(h, k4, k4, PureEndowmentStates);
			saxpy(1 / 6.0f, k4, v, tmp, PureEndowmentStates);
			saxpy(2 / 6.0f, k3, tmp, tmp, PureEndowmentStates);
			saxpy(2 / 6.0f, k2, tmp, tmp, PureEndowmentStates);
			saxpy(1 / 6.0f, k1, tmp, v, PureEndowmentStates);
		}
		floatArrayCopy(v, &result[offset + (y*PureEndowmentStates) - PureEndowmentStates], PureEndowmentStates);
	}
}

__device__ void DeferredTemporaryLifeAnnuity(int a, int b, int steps, float* Va, float* result){
	unsigned long offset = getOffset(threadIdx.x, blockIdx.x, blockDim.x, DeferredTemporaryLifeAnnuityStates, a);

	float h = -1.0f / steps;
	float k1[DeferredTemporaryLifeAnnuityStates], k2[DeferredTemporaryLifeAnnuityStates], k3[DeferredTemporaryLifeAnnuityStates], k4[DeferredTemporaryLifeAnnuityStates], tmp[DeferredTemporaryLifeAnnuityStates], v[DeferredTemporaryLifeAnnuityStates];

	floatArrayCopy(Va, &result[offset + a * DeferredTemporaryLifeAnnuityStates], DeferredTemporaryLifeAnnuityStates);

	for(int y = a; y > b; y--) {
		dtla_bj_ii(y, v);
		saxpy(1.0, v, &result[offset + y * DeferredTemporaryLifeAnnuityStates], v, DeferredTemporaryLifeAnnuityStates);
		for(int s = 0; s < steps; s++) { 	// Integrate backwards over [y,y-1]
			float t = y - s / (float)steps;
			// Hack: Fake limit from left
			dtla_dV(s == 0 ? t - 1e-5 : t, v, k1);
			sax(h, k1, k1, DeferredTemporaryLifeAnnuityStates);
			saxpy(0.5, k1, v, tmp, DeferredTemporaryLifeAnnuityStates);
			dtla_dV(t + h / 2, tmp, k2);
			sax(h, k2, k2, DeferredTemporaryLifeAnnuityStates);
			saxpy(0.5, k2, v, tmp, DeferredTemporaryLifeAnnuityStates);
			dtla_dV(t + h / 2, tmp, k3);
			sax(h, k3, k3, DeferredTemporaryLifeAnnuityStates);
			saxpy(1, k3, v, tmp, DeferredTemporaryLifeAnnuityStates);
			// Hack: Fake limit from right
			dtla_dV(s == steps - 1 ? t + h + 1e-5 : t + h, tmp, k4);
			sax(h, k4, k4, DeferredTemporaryLifeAnnuityStates);
			saxpy(1 / 6.0f, k4, v, tmp, DeferredTemporaryLifeAnnuityStates);
			saxpy(2 / 6.0f, k3, tmp, tmp, DeferredTemporaryLifeAnnuityStates);
			saxpy(2 / 6.0f, k2, tmp, tmp, DeferredTemporaryLifeAnnuityStates);
			saxpy(1 / 6.0f, k1, tmp, v, DeferredTemporaryLifeAnnuityStates);
		}
		floatArrayCopy(v, &result[offset + (y*DeferredTemporaryLifeAnnuityStates) - DeferredTemporaryLifeAnnuityStates], DeferredTemporaryLifeAnnuityStates);
	}
}

__device__ void TemporaryLifeAnnuityPremium(int a, int b, int steps, float* Va, float* result){
	unsigned long offset = getOffset(threadIdx.x, blockIdx.x, blockDim.x, TemporaryLifeAnnuityPremiumStates, a);

	float h = -1.0f / steps;
	float k1[TemporaryLifeAnnuityPremiumStates], k2[TemporaryLifeAnnuityPremiumStates], k3[TemporaryLifeAnnuityPremiumStates], k4[TemporaryLifeAnnuityPremiumStates], tmp[TemporaryLifeAnnuityPremiumStates], v[TemporaryLifeAnnuityPremiumStates];

	floatArrayCopy(Va, &result[offset + a * TemporaryLifeAnnuityPremiumStates], TemporaryLifeAnnuityPremiumStates);

	for(int y = a; y > b; y--) {
		tlap_bj_ii(y, v);
		saxpy(1.0, v, &result[offset + y * TemporaryLifeAnnuityPremiumStates], v, TemporaryLifeAnnuityPremiumStates);
		for(int s = 0; s < steps; s++) { 	// Integrate backwards over [y,y-1]
			float t = y - s / (float)steps;
			// Hack: Fake limit from left
			tlap_dV(s == 0 ? t - 1e-5 : t, v, k1);
			sax(h, k1, k1, TemporaryLifeAnnuityPremiumStates);
			saxpy(0.5, k1, v, tmp, TemporaryLifeAnnuityPremiumStates);
			tlap_dV(t + h / 2, tmp, k2);
			sax(h, k2, k2, TemporaryLifeAnnuityPremiumStates);
			saxpy(0.5, k2, v, tmp, TemporaryLifeAnnuityPremiumStates);
			tlap_dV(t + h / 2, tmp, k3);
			sax(h, k3, k3, TemporaryLifeAnnuityPremiumStates);
			saxpy(1, k3, v, tmp, TemporaryLifeAnnuityPremiumStates);
			// Hack: Fake limit from right
			tlap_dV(s == steps - 1 ? t + h + 1e-5 : t + h, tmp, k4);
			sax(h, k4, k4, TemporaryLifeAnnuityPremiumStates);
			saxpy(1 / 6.0f, k4, v, tmp, TemporaryLifeAnnuityPremiumStates);
			saxpy(2 / 6.0f, k3, tmp, tmp, TemporaryLifeAnnuityPremiumStates);
			saxpy(2 / 6.0f, k2, tmp, tmp, TemporaryLifeAnnuityPremiumStates);
			saxpy(1 / 6.0f, k1, tmp, v, TemporaryLifeAnnuityPremiumStates);
		}
		floatArrayCopy(v, &result[offset + (y*TemporaryLifeAnnuityPremiumStates) - TemporaryLifeAnnuityPremiumStates], TemporaryLifeAnnuityPremiumStates);
	}
}

__device__ void TermInsurance(int a, int b, int steps, float* Va, float* result){
	unsigned long offset = getOffset(threadIdx.x, blockIdx.x, blockDim.x, TermInsuranceStates, a);

	float h = -1.0f / steps;
	float k1[TermInsuranceStates], k2[TermInsuranceStates], k3[TermInsuranceStates], k4[TermInsuranceStates], tmp[TermInsuranceStates], v[TermInsuranceStates];

	floatArrayCopy(Va, &result[offset + a * TermInsuranceStates], TermInsuranceStates);

	for(int y = a; y > b; y--) {
		ti_bj_ii(y, v);
		saxpy(1.0, v, &result[offset + y * TermInsuranceStates], v, TermInsuranceStates);
		for(int s = 0; s < steps; s++) { 	// Integrate backwards over [y,y-1]
			float t = y - s / (float)steps;
			// Hack: Fake limit from left
			ti_dV(s == 0 ? t - 1e-5 : t, v, k1);
			sax(h, k1, k1, TermInsuranceStates);
			saxpy(0.5, k1, v, tmp, TermInsuranceStates);
			ti_dV(t + h / 2, tmp, k2);
			sax(h, k2, k2, TermInsuranceStates);
			saxpy(0.5, k2, v, tmp, TermInsuranceStates);
			ti_dV(t + h / 2, tmp, k3);
			sax(h, k3, k3, TermInsuranceStates);
			saxpy(1, k3, v, tmp, TermInsuranceStates);
			// Hack: Fake limit from right
			ti_dV(s == steps - 1 ? t + h + 1e-5 : t + h, tmp, k4);
			sax(h, k4, k4, TermInsuranceStates);
			saxpy(1 / 6.0f, k4, v, tmp, TermInsuranceStates);
			saxpy(2 / 6.0f, k3, tmp, tmp, TermInsuranceStates);
			saxpy(2 / 6.0f, k2, tmp, tmp, TermInsuranceStates);
			saxpy(1 / 6.0f, k1, tmp, v, TermInsuranceStates);
		}
		floatArrayCopy(v, &result[offset + (y*TermInsuranceStates) - TermInsuranceStates], TermInsuranceStates);
	}
}

__device__ void DisabilityAnnuity(int a, int b, int steps, float* Va, float* result){
	unsigned long offset = getOffset(threadIdx.x, blockIdx.x, blockDim.x, DisabilityAnnuityStates, a);

	float h = -1.0f / steps;
	float k1[DisabilityAnnuityStates], k2[DisabilityAnnuityStates], k3[DisabilityAnnuityStates], k4[DisabilityAnnuityStates], tmp[DisabilityAnnuityStates], v[DisabilityAnnuityStates];

	floatArrayCopy(Va, &result[offset + a * DisabilityAnnuityStates], DisabilityAnnuityStates);

	for(int y = a; y > b; y--) {
		da_bj_ii(y, v);
		saxpy(1.0, v, &result[offset + y * DisabilityAnnuityStates], v, DisabilityAnnuityStates);
		for(int s = 0; s < steps; s++) { 	// Integrate backwards over [y,y-1]
			float t = y - s / (float)steps;
			// Hack: Fake limit from left
			da_dV(s == 0 ? t - 1e-5 : t, v, k1);
			sax(h, k1, k1, DisabilityAnnuityStates);
			saxpy(0.5, k1, v, tmp, DisabilityAnnuityStates);
			da_dV(t + h / 2, tmp, k2);
			sax(h, k2, k2, DisabilityAnnuityStates);
			saxpy(0.5, k2, v, tmp, DisabilityAnnuityStates);
			da_dV(t + h / 2, tmp, k3);
			sax(h, k3, k3, DisabilityAnnuityStates);
			saxpy(1, k3, v, tmp, DisabilityAnnuityStates);
			// Hack: Fake limit from right
			da_dV(s == steps - 1 ? t + h + 1e-5 : t + h, tmp, k4);
			sax(h, k4, k4, DisabilityAnnuityStates);
			saxpy(1 / 6.0f, k4, v, tmp, DisabilityAnnuityStates);
			saxpy(2 / 6.0f, k3, tmp, tmp, DisabilityAnnuityStates);
			saxpy(2 / 6.0f, k2, tmp, tmp, DisabilityAnnuityStates);
			saxpy(1 / 6.0f, k1, tmp, v, DisabilityAnnuityStates);
		}
		floatArrayCopy(v, &result[offset + (y*DisabilityAnnuityStates) - DisabilityAnnuityStates], DisabilityAnnuityStates);
	}
}

__device__ void DisabilityTermInsurance(int a, int b, int steps, float* Va, float* result){
	unsigned long offset = getOffset(threadIdx.x, blockIdx.x, blockDim.x, DisabilityTermInsuranceStates, a);

	float h = -1.0f / steps;
	float k1[DisabilityTermInsuranceStates], k2[DisabilityTermInsuranceStates], k3[DisabilityTermInsuranceStates], k4[DisabilityTermInsuranceStates], tmp[DisabilityTermInsuranceStates], v[DisabilityTermInsuranceStates];

	floatArrayCopy(Va, &result[offset + a * DisabilityTermInsuranceStates], DisabilityTermInsuranceStates);

	for(int y = a; y > b; y--) {
		dti_bj_ii(y, v);
		saxpy(1.0, v, &result[offset + y * DisabilityTermInsuranceStates], v, DisabilityTermInsuranceStates);
		for(int s = 0; s < steps; s++) { 	// Integrate backwards over [y,y-1]
			float t = y - s / (float)steps;
			// Hack: Fake limit from left
			dti_dV(s == 0 ? t - 1e-5 : t, v, k1);
			sax(h, k1, k1, DisabilityTermInsuranceStates);
			saxpy(0.5, k1, v, tmp, DisabilityTermInsuranceStates);
			dti_dV(t + h / 2, tmp, k2);
			sax(h, k2, k2, DisabilityTermInsuranceStates);
			saxpy(0.5, k2, v, tmp, DisabilityTermInsuranceStates);
			dti_dV(t + h / 2, tmp, k3);
			sax(h, k3, k3, DisabilityTermInsuranceStates);
			saxpy(1, k3, v, tmp, DisabilityTermInsuranceStates);
			// Hack: Fake limit from right
			dti_dV(s == steps - 1 ? t + h + 1e-5 : t + h, tmp, k4);
			sax(h, k4, k4, DisabilityTermInsuranceStates);
			saxpy(1 / 6.0f, k4, v, tmp, DisabilityTermInsuranceStates);
			saxpy(2 / 6.0f, k3, tmp, tmp, DisabilityTermInsuranceStates);
			saxpy(2 / 6.0f, k2, tmp, tmp, DisabilityTermInsuranceStates);
			saxpy(1 / 6.0f, k1, tmp, v, DisabilityTermInsuranceStates);
		}
		floatArrayCopy(v, &result[offset + (y*DisabilityTermInsuranceStates) - DisabilityTermInsuranceStates], DisabilityTermInsuranceStates);
	}
}

__global__ void RunPlans(float *peVa, float *peResult, 
						 float *dtlaVa, float *dtlaResult,
						 float *tlapVa, float *tlapResult, 
						 float *tiVa, float *tiResult, 
						 float *daVa, float *daResult, 
						 float *dtiVa, float *dtiResult){
	int steps = 100;

	PureEndowment(40, 0, steps, peVa, peResult);
	DeferredTemporaryLifeAnnuity(50, 0, steps, dtlaVa, dtlaResult);
	TemporaryLifeAnnuityPremium(50, 0, steps, tlapVa, tlapResult);
	TermInsurance(50, 0, steps, tiVa, tiResult);
	DisabilityAnnuity(50, 0, steps, daVa, daResult);
	DisabilityTermInsurance(50, 0, steps, dtiVa, dtiResult);
}

/*
function pointers not allowed on sm_1.x
malloc not allowed on device

typedef void (*dV)(float t, float *v, float *res);
typedef void (*bj_ii)(float t, float *res);

__global__ void RK4_n(dV dV, bj_ii bj_ii, int a, int b, int steps, float* Va, float* result){
unsigned long offset = getOffset(threadIdx.x, blockIdx.x, blockDim.x, PureEndowmentStates, a); // (threadIdx.x + blockIdx.x * blockDim.x) * n * (a + 1)

float h = -1.0f / steps;
float k1[PureEndowmentStates], k2[PureEndowmentStates], k3[PureEndowmentStates], k4[PureEndowmentStates], tmp[PureEndowmentStates], v[PureEndowmentStates];

floatArrayCopy(Va, &result[offset + a * PureEndowmentStates], PureEndowmentStates); // result[a - b] = Va;

for(int y = a; y > b; y--) {
bj_ii(y, v);
saxpy(1.0, v, &result[offset + y * PureEndowmentStates], v, PureEndowmentStates);
for(int s = 0; s < steps; s++) { 	// Integrate backwards over [y,y-1]
float t = y - s / (float)steps;
// Hack: Fake limit from left
dV(s == 0 ? t - 1e-14 : t, v, k1);
sax(h, k1, k1, PureEndowmentStates);
saxpy(0.5, k1, v, tmp, PureEndowmentStates);
dV(t + h / 2, tmp, k2);
sax(h, k2, k2, PureEndowmentStates);
saxpy(0.5, k2, v, tmp, PureEndowmentStates);
dV(t + h / 2, tmp, k3);
sax(h, k3, k3, PureEndowmentStates);
saxpy(1, k3, v, tmp, PureEndowmentStates);
// Hack: Fake limit from right
dV(s == steps - 1 ? t + h + 1e-14 : t + h, tmp, k4);
sax(h, k4, k4, PureEndowmentStates);
saxpy(1 / 6.0f, k4, v, tmp, PureEndowmentStates);
saxpy(2 / 6.0f, k3, tmp, tmp, PureEndowmentStates);
saxpy(2 / 6.0f, k2, tmp, tmp, PureEndowmentStates);
saxpy(1 / 6.0f, k1, tmp, v, PureEndowmentStates);
}
floatArrayCopy(v, &result[offset + (y*PureEndowmentStates) - PureEndowmentStates], PureEndowmentStates); //Array.Copy(v, result[y - 1 - b], v.Length);
}
}
*/