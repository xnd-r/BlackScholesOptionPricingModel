#include "NumSolution.h"

float NumSolution::EulMarStep(float S, float dt, float dw, float r, float sig) {
	return S + S * (r * dt + sig * dw);
}

float NumSolution::MilsteinStep(float S, float dt, float dw, float r, float sig) {
	return S + S * (r * dt + sig * dw) + 0.5f * S * sig * sig * (dw * dw - dt);
}

float NumSolution::RK1Step(float S, float dt, float dw, float r, float sig) {
	return S * (1.0f + r * dt + sig * dw + 0.5f * sig * (r * dt + sig * sqrtf(dt)) * (dw * dw - dt) / sqrtf(dt));
	// Ñonverges with nsamples == 5000
}

float NumSolution::BurragePlatenStep(float S, float dt, float dw, float dz, float r, float sig) {
	return S * (1.f + r * dt + sig * dw) + .5f * S * sig * sig * (dw * dw - dt) +
		r * sig * S * dz + .5f * r * r * S * dt * dt +
		r * sig * S * (dw * dt - dz) + .5f * sig * sig * sig * S *
		(1.f / 3.f * dw * dw - dt) * dw;
	// Ñonverges with nsamples == 135
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


	for (int j = 2; j <= nsteps * 2; j += 2)
	{
		buffer[j] = buffer[j - 2] + dw[j - 2];
		buffer[j + 1] = buffer[j - 1] + dw[j - 1];
	}

	delete[] dw;
}

bool NumSolution::IsConvergence(int StepIndex, VSLStreamStatePtr stream, int npaths, int nsteps, float pS0, float pR, float pSig, float time, float *Error)
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
			scale = scale / 2;
		}
	}

	for (int j = 0; j < 8; j++)
		Error[j] = Error[j] / npaths;
	delete[] wtraject;
	return Error[7] < 0.01 ? true : false;
}

float NumSolution::stockPricesIntegrator(VSLStreamStatePtr stream, float* wtraject, int StepIndex, int nsteps, float pS0, float pR, float pSig, float time)
{
	float GlobalStockPrice = 0.0f;
	float dt = time / nsteps;
	float t = 0.0f;
	float stockPrice;
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
	return stockPrice;
}

float NumSolution::SimulateStockPrices(int StepIndex, int npaths, int nsteps, float pS0, float pR, float pSig, float time)
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

float NumSolution::SimulateStockPricesVol(int StepIndex, int npaths, int nsteps, float pS0, float* pR, float* pSig, float time)
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

float NumSolution::getMCPrice(int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0) {
	VSLStreamStatePtr stream;
	float* wtraject;
	if (StepIndex == 3) {
		wtraject = new float[(nsteps + 1) * 2];
	}
	else {
		wtraject = new float[nsteps + 1];
	}

	if (indexGen == 0) {
		const unsigned int _seed[2] = { seed, seed };
		vslNewStreamEx(&stream, VSL_BRNG_MCG59, 2, _seed);
	}
	else if (indexGen == 1) {
		vslNewStream(&stream, VSL_BRNG_SOBOL, 1);
	}

	float payoff, sum = .0f, stockPrice;
	for (unsigned int portion = 0; portion < N; portion++) {
			stockPrice = stockPricesIntegrator(stream, wtraject, StepIndex, nsteps, pS0, R, SIG, Time);
			payoff = stockPrice - K;
			if (payoff > 0.f)
				sum += payoff;
	}
	sum = sum / N * expf(-R * Time);
	vslDeleteStream(&stream);
	delete[] wtraject;
	return sum;
}

float NumSolution::getMCPricePar(int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0) {
	float* wtraject;
	if (StepIndex == 3) {
		wtraject = new float[(nsteps + 1) * 2];
	}
	else {
		wtraject = new float[nsteps + 1];
	}
	float sum = .0f, stockPrice;
#pragma omp parallel private(stockPrice)
	{
		VSLStreamStatePtr stream;
		int count = omp_get_num_threads();
		int num = omp_get_thread_num();

		if (indexGen == 0) {
			const unsigned int _seed[2] = { seed, seed };
			vslNewStreamEx(&stream, VSL_BRNG_MCG59, 2, _seed);
		}
		else if (indexGen == 1) {
			vslNewStream(&stream, VSL_BRNG_SOBOL, 1);
		}
		vslSkipAheadStream(stream, N / count * num);
		float payoff;
#pragma omp for private(payoff), reduction(+:sum)
		for (unsigned int portion = 0; portion < N; portion++) {
			stockPrice = stockPricesIntegrator(stream, wtraject, StepIndex, nsteps, pS0, R, SIG, Time);
			payoff = stockPrice - K;
			if (payoff > 0.f)
				sum += payoff;
		}
		vslDeleteStream(&stream);
	}
	sum = sum / N * expf(-R * Time);

	delete[] wtraject;
	return sum;
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
		fprintf(f, "%stockPrice\n", tmp_string.c_str());
	}
	fclose(f);
}

void NumSolution::getErrors(int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0, int sampleStep) {
	time_t timestamp;
	time(&timestamp);
	std::string date = asctime(localtime(&timestamp));
	date.pop_back();
	std::string fileName;
	date.append(std::to_string(N));
	date.append("_Samples_");
	date.append(std::to_string(sampleStep));
	date.append("_Step_");
	fileName = "_Errors.csv";
	date.append(fileName);


	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}


	FILE *f = fopen(date.c_str(), "w");

	fprintf(f, "%stockPrice;%stockPrice;%stockPrice;%stockPrice\n;", "Euler", "Mils", "RK1", "Platen");
	for (int i = 0; i < N; i += sampleStep) {
		fprintf(f, "%lf;%lf;%lf;%lf\n", 
			(getMCPrice(0, nsteps, indexGen, i, seed, K, R, Time, SIG, pS0) - 20.9244),
			(getMCPrice(1, nsteps, indexGen, i, seed, K, R, Time, SIG, pS0) - 20.9244),
			(getMCPrice(2, nsteps, indexGen, i, seed, K, R, Time, SIG, pS0) - 20.9244), 
			(getMCPrice(3, nsteps, indexGen, i, seed, K, R, Time, SIG, pS0) - 20.9244));
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


float NumSolution::Integrand(float z) {
	float payoff;	
	int scale = 2000;
	stockPrice = pS0 * expf(tmp1 + tmp2 * z);
	payoff = stockPrice - K;
	payoff > 0.0f ? payoff *= expf(-z * z / 2.0f) : payoff = 0.0f;
	return payoff;
}

float NumSolution::GetRPrice(float a, float b) {
	int scale = 2000;
	float h = (b - a) / scale;
	start = clock();
	for (int i = 0; i < scale; ++i)
		sum += Integrand(a + h * (i + 0.5f));
	sum = sum * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * h;

	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;

	return sum;
}

float NumSolution::GetTPrice(float a, float b) {
	int scale = 2000;
	float h = (b - a) / scale;
	start = clock();

	sum += Integrand(a) * 0.5f * h + Integrand(b) * 0.5f * (b - a);
	for (int i = 1; i < scale - 1; ++i)
		sum += Integrand(a + h * i);
	sum = sum * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * h;

	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;

	return sum;
}

float NumSolution::GetSPrice(float a, float b) {
	int scale = 2000;
	float h = (b - a) / scale;
	float sum2 = 0.0f;

	start = clock();
	float sum4 = Integrand(a + h);
	sum += Integrand(a) + Integrand(b);

	for (int i = 1; i < scale - 2; i += 2) {
		sum4 += Integrand(a + (i + 1)*h);
		sum2 += Integrand(a + i * h);
	}

	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;
	sum = (sum + 4 * sum4 + 2 * sum2) * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * (h / 3.0f);

	return sum;
}

float NumSolution::Get3_8Price(float a, float b) {
	int scale = 2000;
	float h = (b - a) / scale;
	float sum2 = 0.0f;
	float sum3 = Integrand(a + h);

	start = clock();
	sum += Integrand(a) + Integrand(b);

	for (int i = 1; i < scale - 2; i += 3) {
		sum3 += Integrand(a + i * h);
	}

	for (int i = 2; i < scale - 1; i += 3) {
		sum3 += Integrand(a + i * h);
	}

	for (int i = 3; i < scale; i += 3) {
		sum2 += Integrand(a + i * h);
	}

	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;
	sum = (sum + 3 * sum3 + 2 * sum2) * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * (3.0f * h / 8.0f);

	return sum;
}

void NumSolution::SetS(int amo) {
	int scale = 2000;
	s_array = new float[amo];
	exp_array = new float[amo];
#if defined(__INTEL_COMPILER) 
#pragma simd
#pragma vector always	
#endif
	int j = 0;
#pragma omp parallel for private(j) 
	for (j = 0; j < amo; ++j) {
		s_array[j] = pS0 * expf(tmp1 + tmp2 * (a + h * ((float)j + 0.5f)));
		exp_array[j] = expf(-(a + h * ((float)j + 0.5f)) * (a + h * ((float)j + 0.5f)) / 2.0f);
	}
}

float NumSolution::GetRPrice(float a, float b, int NumThreads/*, float* s_array, float* expf_array*/) {
	float sum = 0.0f;
	int scale = 2000;
	float h = (b - a) / scale;
	float payoff;
	int i, ind;
	start = clock();
	omp_set_num_threads(NumThreads);
#if defined(__INTEL_COMPILER) 
#pragma simd
#pragma vector always	
#endif
#pragma omp parallel for private(i, payoff, ind) reduction(+:sum)
	for (ind = 0; ind < N; ind++) {
		for (i = 0; i < scale; ++i) {
			payoff = s_array[i] - K;
			payoff > 0.0f ? payoff *= exp_array[i] : payoff = 0.0f;
			sum = sum + payoff;
		}
	}
	sum = (sum * expf(-R * TIME) * h) / sqrtf(2.0f * M_PIF);
	sum /= N;
	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;

	return sum;
}

float NumSolution::GetTPrice(float a, float b, int NumThreads) {
	float sum;
	int scale = 2000;
	float h = (b - a) / scale;
	float payoff;
	int i, ind;
	start = clock();
	omp_set_num_threads(NumThreads);
#if defined(__INTEL_COMPILER) 
#pragma simd
#pragma vector always	
#endif

	sum = (Integrand(a) + Integrand(b)) * 0.5f * h * (b - a);
#pragma omp parallel for private(i, payoff, ind) reduction(+:sum)
	for (ind = 0; ind < N; ind++) {
		for (int i = 1; i < scale - 1; ++i) {
			payoff = s_array[i] - K;
			payoff > 0.0f ? payoff *= exp_array[i] : payoff = 0.0f;
			sum = sum + payoff;
		}
	}
	sum = sum * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * h;
	sum /= N;
	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;

	return sum;
}

float NumSolution::GetSPrice(float a, float b, int NumThreads) {
	float sum;
	int scale = 2000;
	float h = (b - a) / scale;
	float payoff;
	int i, ind;
	float sum2 = 0.0f;
	start = clock();
	omp_set_num_threads(NumThreads);
#if defined(__INTEL_COMPILER) 
#pragma simd
#pragma vector always	
#endif
	float sum4 = Integrand(a + h);
	sum = Integrand(a) + Integrand(b);

#pragma omp parallel for private(i, payoff, ind) reduction(+:sum2, sum4)
	for (ind = 0; ind < N; ind++) {
		for (int i = 1; i < scale - 2; i += 2) {
			payoff = s_array[i] - K;
			payoff > 0.0f ? payoff *= exp_array[i] : payoff = 0.0f;
			sum2 = sum2 + payoff;

			payoff = s_array[i + 1] - K;
			payoff > 0.0f ? payoff *= exp_array[i + 1] : payoff = 0.0f;
			sum4 = sum4 + payoff;
		}
	}
	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;
	sum = (sum + 4 * sum4 + 2 * sum2) * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * (h / 3.0f);
	sum /= N;
	return sum;
}


float NumSolution::Get3_8Price(float a, float b, int NumThreads) {
	float sum;
	int scale = 2000;
	float h = (b - a) / scale;
	float payoff;
	int i, ind;
	start = clock();
	omp_set_num_threads(NumThreads);
#if defined(__INTEL_COMPILER) 
#pragma simd
#pragma vector always	
#endif
	float sum2 = 0.0f;
	float sum3 = Integrand(a + h);
	sum += Integrand(a) + Integrand(b);

#pragma omp parallel for private(i, payoff, ind) reduction(+:sum2, sum3)
	for (ind = 0; ind < N; ind++) {
		for (int i = 1; i < scale - 2; i += 3) {
			payoff = s_array[i] - K;
			payoff > 0.0f ? payoff *= exp_array[i] : payoff = 0.0f;
			sum3 = sum3 + payoff;
		}

		for (int i = 2; i < scale - 1; i += 3) {
			payoff = s_array[i] - K;
			payoff > 0.0f ? payoff *= exp_array[i] : payoff = 0.0f;
			sum3 = sum3 + payoff;
		}

		for (int i = 3; i < scale; i += 3) {
			payoff = s_array[i] - K;
			payoff > 0.0f ? payoff *= exp_array[i] : payoff = 0.0f;
			sum2 = sum2 + payoff;
		}
	}
	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;
	sum = (sum + 3 * sum3 + 2 * sum2) * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * (3.0f * h / 8.0f);
	sum /= N;
	return sum;
}

void NumSolution::Execute() {
	int tmp = 1;
	double t1, t2;
	float tmp_price;
	int scale = 2000;
	SetS(scale);
	for (int k = 1; k < 3; ++k) {
		for (int j = 0; j < 7; ++j) {
			t1 = omp_get_wtime();
			//for (int j = 0; j < N; ++j)
			tmp_price = Get3_8Price(-5.15f, 6.0f, tmp);
			t2 = omp_get_wtime();
			//std::cout << j << " Launch Done " << std::endl;
			Times.push_back(t2 - t1);
			Prices.push_back(tmp_price);
			std::cout << j << " Writing Done " << std::endl;
		}
		std::sort(Times.begin(), Times.end());
		for (std::vector<double>::const_iterator it = Times.begin(); it != Times.end(); ++it)
			std::cout << *it << std::endl;
		std::cout << " ####################### " << std::endl;
		for (std::vector<double>::const_iterator it = Prices.begin(); it != Prices.end(); ++it)
			std::cout << *it << std::endl;
		std::cout << tmp << " thread(stockPrice) done" << std::endl;
		tmp *= 2;
		Times.clear();
		Prices.clear();
	}
	delete[] s_array;
	delete[] exp_array;

}