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
	for (int i = 0; i < scale - 1; ++i) 
		sum += Integrand(a + h * (i + 0.5f));
	sum = sum * expf(-R * TIME) / (SIG * sqrtf(2.0f * M_PIF));

	return sum;
}