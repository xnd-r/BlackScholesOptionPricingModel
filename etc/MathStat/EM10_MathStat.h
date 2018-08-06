#ifndef ____STATISTICS____
#define ____STATISTICS____

#define _EPSILON_	5E-03f // calculations accuracy

#include <vector>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <string>

#include "../../modules/StockPriece/include/S10_Wiener.h"

class RVCharacteristics
{
public:
	float* wd_sorted;
	int len;

	float MeanCh;
	float VarCh;
	float MedCh;

	float MeanTh;
	float VarTh;
	float MedTh;

	float VarDiff;

public:
	RVCharacteristics(float* wd, int _len, float _h);
	// TODO: add gist, output(everything), graph, Chi-squared...  

	float GetMean();
	float GetVariance();
	float GetMediane();
	bool IsExpValCorrect(float _eps);
	bool IsVArianceCorrect(float var, float _eps);
	void QuickSort(float* array, int first, int last);
	void WriteToCsv(float eps);
	void Execute();
};


#endif // !____STATISTICS____