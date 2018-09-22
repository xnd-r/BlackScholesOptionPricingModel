#include "O22_QuadFormula.h"

// TODO: O21_MonteCarlo.cpp(50): warning C4101: 'i, payoff': unreferenced local variables, but 'ind' is OK

float QuadratureFormula::Integrand(float z) {
	float payoff;
	s = S0 * expf(tmp1 + tmp2 * z);
	payoff = s - K;
	payoff > 0.0f ? payoff *= expf(-z * z / 2.0f) : payoff = 0.0f;
	return payoff;
}

float QuadratureFormula::GetRPrice(float a, float b) {
	h = (b - a) / scale;
	start = clock();
	for (int i = 0; i < scale; ++i)
		sum += Integrand(a + h * (i + 0.5f));
	sum = sum * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * h;

	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;

	return sum;
}

float QuadratureFormula::GetTPrice(float a, float b) {
	h = (b - a) / scale;
	start = clock();

	sum += Integrand(a) * 0.5f * h + Integrand(b) * 0.5f * (b - a);
	for (int i = 1; i < scale - 1; ++i)
		sum += Integrand(a + h * i);
	sum = sum * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * h;

	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;

	return sum;
}

float QuadratureFormula::GetSPrice(float a, float b) {
	h = (b - a) / scale;
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

float QuadratureFormula::Get3_8Price(float a, float b) {
	h = (b - a) / scale;
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

void QuadratureFormula::SetS(int amo) {
	s_array = new float[amo];
	exp_array = new float[amo];
#if defined(__INTEL_COMPILER) 
	#pragma simd
	#pragma vector always	
#endif
	int j = 0;
#pragma omp parallel for private(j) 
	for (j = 0; j < amo; ++j) {
		s_array[j] = S0 * expf(tmp1 + tmp2 * (a + h * ((float)j + 0.5f)));
		exp_array[j] = expf(-(a + h * ((float)j + 0.5f)) * (a + h * ((float)j + 0.5f)) / 2.0f);
	}
}

float QuadratureFormula::GetRPrice(float a, float b, int NumThreads/*, float* s_array, float* expf_array*/) {
	float sum = 0.0f;
	h = (b - a) / scale;
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
	sum = (sum * expf(-R * TIME) * h ) / sqrtf(2.0f * M_PIF);
	sum /= N;
	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;

	return sum;
}

float QuadratureFormula::GetTPrice(float a, float b, int NumThreads) {
	float sum;
	h = (b - a) / scale;
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

float QuadratureFormula::GetSPrice(float a, float b, int NumThreads) {
	float sum;
	h = (b - a) / scale;
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


float QuadratureFormula::Get3_8Price(float a, float b, int NumThreads) {
	float sum;
	h = (b - a) / scale;
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

void QuadratureFormula::Execute() {
	int tmp = 1;
	double t1, t2;
	float tmp_price;
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
		std::cout << tmp << " thread(s) done" << std::endl;
		tmp *= 2;
		Times.clear();
		Prices.clear();
	}
	delete[] s_array;
	delete[] exp_array;

}