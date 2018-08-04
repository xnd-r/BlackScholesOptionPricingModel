#ifndef ____STATISTICS____
#define ____STATISTICS____

#define _EPSILON_	5E-03 // calculations accuracy

#include <vector>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <string>

#include "../../modules/StockPriece/include/S10_Wiener.h"

class RVCharacteristics
{
public:
	double* wd_sorted;
	int len;

	double MeanCh;
	double VarCh;
	double MedCh;

	double MeanTh;
	double VarTh;
	double MedTh;

	double VarDiff;

public:
	RVCharacteristics(double* wd, int _len, double _h);
	// TODO: add gist, output(everything), graph, Chi-squared...  

	double GetMean();
	double GetVariance();
	double GetMediane();
	bool IsExpValCorrect(double _eps);
	bool IsVArianceCorrect(double var, double _eps);
	void QuickSort(double* array, int first, int last);
	void WriteToCsv(double eps);
	void Execute();
};


#endif // !____STATISTICS____