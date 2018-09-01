//pattern for class

#ifndef ____QUADRATURE_FORMULA____
#define ____QUADRATURE_FORMULA____
#include "../NumSolution/O20_NumSolution.h"

#define M_PIF 3.14159265358979323846f

class QuadratureFormula : public NumSolutionOption {
	int scale = 2048; // amount of segments in sum of Darbu
	float a = -5.15f;
	float b = 6.0f;
	float h = (b - a) / scale;
	public:

	typedef float(QuadratureFormula::*Method)(float, float);

	QuadratureFormula::Method MethodArray[3] = { &QuadratureFormula::GetRPrice, &QuadratureFormula::GetTPrice, &QuadratureFormula::GetSPrice };

	float Integrand(float x);
		
	float GetRPrice	(float a, float b);
	float GetTPrice(float a, float b);
	float GetSPrice(float a, float b);
	float Get3_8Price(float a, float b);
	void SetS(int amo);
	void Execute();

	float* s_array;
	float* exp_array;

	float GetRPrice(float a, float b, int NumThread/*, float* s_array, float* expf_array*/);
	float GetTPrice(float a, float b, int NumThread);
	float GetSPrice(float a, float b, int NumThread);
	float Get3_8Price(float a, float b, int NumThread);

	

};
#endif // !____QUADRATURE_FORMULA____