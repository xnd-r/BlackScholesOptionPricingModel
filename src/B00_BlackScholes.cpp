#include "../include/B00_BlackScholes.h"

float BlackScholes::GetStockPrice(float wiener_diff, float time) {
	return S0 * expf((R - SIG * SIG / 2.0f) * time + SIG * wiener_diff);
}

