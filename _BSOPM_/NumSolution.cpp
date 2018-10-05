#include "NumSolution.h"

float NumSolution::EulMarStep(float S, float dt, float dw, float r, float sig) {
	return S + S * (r * dt + sig * dw);
}

float NumSolution::MilsteinStep(float S, float dt, float dw, float r, float sig) {
	return S + S * (r * dt + sig * dw) + 0.5f * S * sig * sig * (dw * dw - dt);
}

float NumSolution::RK1Step(float S, float dt, float dw, float r, float sig) {
	return S * (1.0f + r * dt + sig * dw + 0.5f * sig * (r * dt + sig * sqrtf(dt)) * (dw * dw - dt) / sqrtf(dt));
}

float NumSolution::BurragePlatenStep(float S, float dt, float dw, float dz, float r, float sig) {
	float res = S * (1.f + r * dt + sig * dw) + .5f * S * sig * sig * (dw * dw - dt) +
		r * sig * S * dz + .5f * r * r * S * dt * dt +
		r * sig * S * (dw * dt - dz) + .5f * sig * sig * sig * S *
		(1.f / 3.f * dw * dw - dt) * dw;
	return res;
}

void NumSolution::wAndZProcesses(VSLStreamStatePtr stream, int nsteps, float time, float *buffer) {
	float *dw = new float[nsteps * 2];

	float dt = time / nsteps;

	// Getting nsteps random vectors, dim == 2
	float mean[2]; mean[0] = 0.0f; mean[1] = 0.0f;
	float hh = dt * sqrtf(dt);
	float cov[3]; cov[0] = sqrtf(dt); cov[1] = 0.5f * hh;
	cov[2] = 1.0f / sqrtf(12.0f) * hh;

	vsRngGaussianMV(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nsteps, dw, 2, VSL_MATRIX_STORAGE_PACKED, mean, cov);
	buffer[0] = 0; buffer[1] = 0;

#if defined(__INTEL_COMPILER) 
#pragma simd
#endif

	for (int j = 2; j <= nsteps * 2; j += 2)
	{
		buffer[j] = buffer[j - 2] + dw[j - 2];
		buffer[j + 1] = buffer[j - 1] + dw[j - 1];
	}

	delete[] dw;
}

bool NumSolution::IsConvergence(int StepIndex, VSLStreamStatePtr stream, int npaths, int nsteps, 
									float pS0, float pR, float pSig, float time, float *Error)
{	
	for (int i = 0; i < 8; i++)
		Error[i] = 0.0f;
	float* wtraject;
	if (StepIndex != 3) {
		wtraject = new float[nsteps + 1];
	}
	else {
		wtraject = new float[(nsteps + 1) * 2];
	}
	float S_An;
	for (int i = 0; i < npaths; i++) {

		if (StepIndex != 3) {
			wienerProcess(stream, nsteps, time, wtraject);
			S_An = getStockPrice(pS0, pR, pSig, wtraject[nsteps], time);
		}
		else {
			wAndZProcesses(stream, nsteps, time, wtraject);
			S_An = getStockPrice(pS0, pR, pSig, wtraject[nsteps * 2], time);
		}
		float S_Num[8]; // array of S(t) for different scaling 
		int scale = 128;

		for (int j = 0; j < 8; j++) {
			S_Num[j] = pS0;
			int numMethodSteps = nsteps / scale;
			float dt = time / numMethodSteps;
			float t = 0.0f;
			int index = 0;
			switch (StepIndex) {
			case 0: // EulMarStep
				for (int k = 0; k < numMethodSteps; k++) {
					t = t + dt;
					index += scale;
					S_Num[j] = EulMarStep(S_Num[j], dt, wtraject[index] - wtraject[index - scale], pR, pSig);
				}
				break;
			case 1: // MilsteinStep
				for (int k = 0; k < numMethodSteps; k++) {
					t = t + dt;
					index += scale;
					S_Num[j] = MilsteinStep(S_Num[j], dt, wtraject[index] - wtraject[index - scale], pR, pSig);
				}
				break;
			case 2: // RK1Step
				for (int k = 0; k < numMethodSteps; k++) {
					t = t + dt;				
					index += scale;
					S_Num[j] = RK1Step(S_Num[j], dt, wtraject[index] - wtraject[index - scale], pR, pSig);
				}
				break;
			case 3: // BurragePlatenStep
				float dw, dz;
				for (int k = 0; k < numMethodSteps; k++) {
					t = t + dt;
					index += scale;
					dw = wtraject[index * 2] - wtraject[(index - scale) * 2];
					dz = wtraject[index * 2 + 1] - wtraject[(index - scale) * 2 + 1];
					S_Num[j] = BurragePlatenStep(S_Num[j], dt, dw, dz, pR, pSig);
				}
				break;
			}
			Error[j] += fabs(S_An - S_Num[j]);
			//printf("Error %d = %lf\n", j + 1, Error[j]);
			scale = scale / 2;
		}
	}

	for (int j = 0; j < 8; j++)
		Error[j] = Error[j] / npaths;
	delete[] wtraject;
	return Error[7] < 0.01 ? true : false;
}

float NumSolution::SimulateStockPrices(int StepIndex, int npaths, int nsteps,
									float pS0, float pR, float pSig, float time)
{
	VSLStreamStatePtr stream = initGen();
	float* wtraject;
	if (StepIndex == 3) {
		wtraject = new float[(nsteps + 1) * 2];
	}
	else {
		wtraject = new float[nsteps + 1];
	}
	float GlobalStockPrice = 0.0f;
	float dt = time / nsteps;
	float t = 0.0f;
	float stockPrice;
	for (int i = 0; i < npaths; i++) {
		stockPrice = pS0;
		if (StepIndex == 3) {
			wAndZProcesses(stream, nsteps, time, wtraject);
		}
		else {
			wienerProcess(stream, nsteps, time, wtraject);
		}


		switch (StepIndex) {
		case 0: // EulMarStep
			for (int k = 0; k < nsteps; k++) {
				t = t + dt;
				stockPrice = EulMarStep(stockPrice, dt, wtraject[k + 1] - wtraject[k], pR, pSig);
			}
			break;
		case 1: // MilsteinStep
			for (int k = 0; k < nsteps; k++) {
				t = t + dt;
				stockPrice = MilsteinStep(stockPrice, dt, wtraject[k + 1] - wtraject[k], pR, pSig);
			}
			break;
		case 2: // RK1Step
			for (int k = 0; k < nsteps; k++) {
				t = t + dt;
				stockPrice = RK1Step(stockPrice, dt, wtraject[k + 1] - wtraject[k], pR, pSig);
			}
			break;
		case 3: // BurragePlatenStep
			float dw, dz;
			for (int k = 0; k < nsteps; k++) {
				t = t + dt;
				dw = wtraject[k * 2 + 2] - wtraject[k * 2];
				dz = wtraject[k * 2 + 3] - wtraject[k * 2 + 1];
				stockPrice = BurragePlatenStep(stockPrice, dt, dw, dz, pR, pSig);
			}
			break;
		}
		GlobalStockPrice += stockPrice;
	}
	delete[] wtraject;
	freeGen(stream);
	return GlobalStockPrice / npaths;
}

float NumSolution::SimulateStockPricesVol(int StepIndex, int npaths, int nsteps,
	float pS0, float* pR, float* pSig, float time)
{
	VSLStreamStatePtr stream = initGen();
	float *wtraject = new float[nsteps + 1];
	float GlobalStockPrice = 0.0f;
	float dt = time / nsteps;
	float t = 0.0f;
	float stockPrice;
	for (int i = 0; i < npaths; i++) {
		stockPrice = pS0;
		wienerProcess(stream, nsteps, time, wtraject);
		switch (StepIndex) {
		case 0: // EulMarStep
			for (int k = 0; k < nsteps; k++) {
				t = t + dt;
				stockPrice = EulMarStep(stockPrice, dt, wtraject[k + 1] - wtraject[k], pR[k], pSig[k]);
			}
			break;
		case 1: // MilsteinStep
			for (int k = 0; k < nsteps; k++) {
				t = t + dt;
				stockPrice = MilsteinStep(stockPrice, dt, wtraject[k + 1] - wtraject[k], pR[k], pSig[k]);
			}
			break;
		case 2: // RK1Step
			for (int k = 0; k < nsteps; k++) {
				t = t + dt;
				stockPrice = RK1Step(stockPrice, dt, wtraject[k + 1] - wtraject[k], pR[k], pSig[k]);
			}
			break;
		case 3: // BurragePlatenStep
			float dw, dz;
			for (int k = 0; k < nsteps; k++) {
				t = t + dt;
				dw = wtraject[k * 2 + 2] - wtraject[k * 2];
				dz = wtraject[k * 2 + 3] - wtraject[k * 2 + 1];
				stockPrice = BurragePlatenStep(stockPrice, dt, dw, dz, pR[k], pSig[k]);
			}
			break;
		}
	GlobalStockPrice += stockPrice;
	}
	delete[] wtraject;
	freeGen(stream);
	return GlobalStockPrice / npaths;
}


void NumSolution::WriteToCsv(float* Errors, int nsteps, int nrows, float Time, int scale, int stepIndex) {
	row_table rt;
	time_t timestamp;
	time(&timestamp);
	std::string date = asctime(localtime(&timestamp));
	date.pop_back();
	std::string fileName;
	switch (stepIndex) {
	case 0:
		fileName = "_EulerMar.csv";
		break;
	case 1:
		fileName = "_Milstein.csv";
		break;
	case 2:
		fileName = "_RK1.csv";
		break;
	case 3:
		fileName = "_BurragePlaten.csv";
		break;

	}
	date.append(fileName);

	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}


	FILE *f = fopen(date.c_str(), "w");
	fprintf(f, "step;e;log(step);log(e);\n");
	for (int i = 0; i < nrows; ++i) {
		rt.setValues(Errors[i], nsteps, static_cast<int>(scale / pow(2, i)), nrows, Time);
		std::string tmp_string = rt.createString();
		for (std::string::iterator it = tmp_string.begin(); it<tmp_string.end(); ++it) {
			std::replace(tmp_string.begin(), tmp_string.end(), '.', ',');
		}
		fprintf(f, "%s\n", tmp_string.c_str());
	}
	fclose(f);
}

void NumSolution::Execute(int StepIndex, int npaths, int nsteps, float pS0, float pR, float pSig, float time) {
	VSLStreamStatePtr stream = initGen();
	float *Error = new float[8];
	IsConvergence(StepIndex, stream, npaths, nsteps, pS0, pR, pSig, time, Error);
	WriteToCsv(Error, nsteps, 8, time, 128, StepIndex);
	for (int i = 0; i < 8; i++)
		printf("Error %d = %lf\n", i + 1, Error[i]);
	freeGen(stream);
	delete[] Error;
}