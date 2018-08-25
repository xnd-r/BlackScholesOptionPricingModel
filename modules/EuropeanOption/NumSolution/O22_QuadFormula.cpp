#include "O22_QuadFormula.h"

float QuadratureFormula::Integrand(float z) {
			float payoff;
			s = S0 * expf(tmp1 + tmp2 * z);
			payoff = s - K;
			if (payoff > 0.0f)
				payoff *= expf(-z * z / 2.0f);
			else payoff = 0.0f;
		
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

	sum += Integrand(a) * 0.5f * h + Integrand(b) * 0.5 * (b - h);
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