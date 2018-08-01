#include "../include/BlackScholes.h"

double BlackSñholes::GetStockPrice(double wiener_diff, double time) {
	return S0 * exp((R - SIG * SIG / 2.0) * time + SIG * wiener_diff);
}

