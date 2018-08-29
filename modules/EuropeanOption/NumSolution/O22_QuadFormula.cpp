#include "O22_QuadFormula.h"

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
	for (int i = 0; i < scale - 1; ++i)
		sum += Integrand(a + h * (i + 0.5f));
	sum = sum * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * h;

	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;

	return sum;
}

float QuadratureFormula::GetTPrice(float a, float b) {
	h = (b - a) / scale;
	start = clock();

	sum += Integrand(a) * 0.5f * h + Integrand(b) * 0.5f * (b - h);
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

	// should we split cycle on two for parallelizing?
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

void QuadratureFormula::SetS(int amo){
	s_array = new float[amo];
	for (int k = 0; k < amo; ++k) {
		s_array[k] = S0 * expf(tmp1 + tmp2 * (a + h * (k + 0.5f)));
	}
}

float QuadratureFormula::GetRPrice(float a, float b, int NumThreads) {
	float sum = 0.0f;
	float s;
	h = (b - a) / scale;
	float payoff;
	int i, k, portion;
	assert(N % bufsize == 0);
	start = clock();
//#pragma omp parallel private(s, portion, i)
//	{
//#pragma omp parallel for private(s, i) reduction(+:sum)
		for (portion = 0; portion < N / bufsize; portion++) {
			for (k = 0; k < bufsize; k++) {
				for (i = 0; i < scale - 1; ++i)
				{
					//s = S0 * expf(tmp1 + tmp2 * (a + h * (i + 0.5f)));
					payoff = s_array[i] - K;
					//std::cout << payoff << " Done " << std::endl;
					if (payoff > 0.0f)
						payoff *= expf(-(a + h * (i + 0.5f)) * (a + h * (i + 0.5f)) / 2.0f);
						else payoff = 0.0f;

					sum += payoff;
				}
			}
		}
//	}
	sum = sum * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * h / N ;

	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;

	return sum;
}

void QuadratureFormula::Execute() {
	int tmp = 4;
	double t1, t2;
	float tmp_price;
	SetS(scale - 1);
	for (int k = 1; k < 2; ++k) {
		omp_set_num_threads(1);
		for (int j = 0; j < 5; ++j) {
			t1 = omp_get_wtime();
			tmp_price = GetRPrice(-5.15f, 6.0f, tmp);
			t2 = omp_get_wtime();
			std::cout << j << " Done " << std::endl;
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

}