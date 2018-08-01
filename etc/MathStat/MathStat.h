#ifndef ____STATISTICS____
#define ____STATISTICS____

#include <vector>
#include <time.h>

class RVCharacteristics
{
	double* wd_sorted;
	// values?

public:
	RVCharacteristics(double* wd, int len);

	// TODO: add mediane, sorting, gist, output(everything), graph, Chi-squared...  

	double GetExpectedVal(double* wd, int len);
	double GetVariance(double* wd, double eval, int len);
	bool IsExpValCorrect(double ev, double _eps);
	bool IsVArianceCorrect(double var, double var_stat, double _eps);
	void QuickSort(double* array, int first, int last);
};


#endif // !____STATISTICS____