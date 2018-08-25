//pattern for class

#ifndef ____QUADRATURE_FORMULA____
#define ____QUADRATURE_FORMULA____
#include "../NumSolution/O20_NumSolution.h"

#define M_PIF 3.14159265358979323846f

class QuadratureFormula : public NumSolutionOption {
	int scale = 1024; // amount of segments in sum of Darbu
	float h; // == (b - a) / scale;
	public:

	typedef float(QuadratureFormula::*Method)(float, float);

	QuadratureFormula::Method MethodArray[3] = { &QuadratureFormula::GetRPrice, &QuadratureFormula::GetTPrice, &QuadratureFormula::GetSPrice };

	float Integrand(float x);
		
	float GetRPrice	(float a, float b);
	float GetTPrice(float a, float b);
	float GetSPrice(float a, float b);
	float Get3_8Price(float a, float b);

	//float GetRPrice(float a, float b, int NumThread) {};
	//float GetTPrice(float a, float b, int NumThread) {};
	//float GetSPrice(float a, float b, int NumThread) {};

	

};
#endif // !____QUADRATURE_FORMULA____