#include "../include/S32_NumMethodsWZ.h"

float NumMethodWZ::BurragePlatenStep(float S, float dt, float dw, float dz) {
	return S * (1.0f + R * dt + SIG * dw) + 0.5f * S * SIG * SIG * (dw * dw - dt) +
		R * SIG * S * dz + 0.5f * R * R * S * dt * dt +
		R * SIG * S * (dw * dt - dz) + 0.5f * SIG * SIG * SIG * S *
		(1.0f / 3.0f * dw * dw - dt) * dw;
}

float NumMethodWZ::Taylor2Step(float S, float dt, float dw, float dz) {
	return S * (1.0f + R * dt + SIG * dw + 0.5f * SIG * SIG * (dw * dw - dt) +
		R * SIG * dz + 0.5f * R * R * dt * dt + R * SIG * (dw * dt - dz));
}

__declspec(noinline) void NumMethodWZ::SimulateWandZProcesses(VSLStreamStatePtr stream, int nSteps, float Time, float *buffer) {
	float *dw = new float[nSteps * 2];

	float h = Time / (float)nSteps;

	// Getting nSteps random vectors, dim == 2
	float mean[2]; mean[0] = 0.0f; mean[1] = 0.0f;
	float hh = h * sqrtf(h);
	float cov[3]; cov[0] = sqrtf(h); cov[1] = 0.5f * hh;
	cov[2] = 1.0f / sqrtf(12.0f) * hh;

	vsRngGaussianMV(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nSteps, dw, 2, VSL_MATRIX_STORAGE_PACKED, mean, cov);

	buffer[0] = 0; buffer[1] = 0;

	#if defined(__INTEL_COMPILER) 
		#pragma simd
		#pragma vector always	
	#endif

	for (int j = 2; j <= nSteps * 2; j += 2)
	{
		buffer[j] = buffer[j - 2] + dw[j - 2];
		buffer[j + 1] = buffer[j - 1] + dw[j - 1];
	}

	delete[] dw;
}

__declspec(noinline) void NumMethodWZ::SimulateStockPrices(Step _step, VSLStreamStatePtr stream, int nPaths, int nSteps, float Time, float *Error) {
	float *w_traject = new float[(nSteps + 1) * 2];

	#if defined(__INTEL_COMPILER) 
		#pragma simd
		#pragma vector always	
	#endif


	for (int i = 0; i < nPaths; i++) {
		SimulateWandZProcesses(stream, nSteps, Time, w_traject);

		float S_An = GetStockPrice(w_traject[nSteps * 2], Time);

		float S_Num[8]; // array of S(t) for different scaling 
		int scale = 128;

		for (int j = 0; j < 8; j++) {
			S_Num[j] = S0;
			int numMethodSteps = nSteps / scale;
			float h = Time / (float)numMethodSteps;
			float t = 0.0f;
			int index = 0;
			for (int k = 0; k < numMethodSteps; k++) {
				t = t + h;
				index = index + scale;
				float dw = w_traject[index * 2] - w_traject[(index - scale) * 2]; 
				float dz = w_traject[index * 2 + 1] - w_traject[(index - scale) * 2 + 1];
				S_Num[j] = (this->*_step)(S_Num[j], h, dw, dz);
			}

			Error[j] = Error[j] + fabs(S_An - S_Num[j]);
			scale = scale / 2;
		}
	}

	for (int j = 0; j < 8; j++)
		Error[j] = Error[j] / nPaths;
	delete[] w_traject;
}

void NumMethodWZ::Execute(Step _step, char* FileName) {
	VSLStreamStatePtr stream = InitGen();
	float *Error = new float[8];
	for (int i = 0; i < 8; i++)
		Error[i] = 0.0f;

	SimulateStockPrices(_step, stream, NPATHS, NSTEPS, TIME, Error);
	WriteToCsv(Error, NSTEPS, 8, TIME, 128, FileName);
	for (int i = 0; i < 8; i++)
		printf("Error %d = %lf\n", i + 1, Error[i]);
	delete[] Error;
	FreeGen(stream);
}