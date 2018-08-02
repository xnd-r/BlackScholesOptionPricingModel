#include "../modules/StockPriece/include/S10_Wiener.h"
#include "../modules/StockPriece/include/S20_AnSimple.h"
#include "../modules/StockPriece/include/S21_AnExtended.h"
#include "../modules/StockPriece/include/S31_NumMethodsW.h"	
#include "../modules/StockPriece/include/S32_NumMethodsWZ.h"	
#include "../etc/MathStat/MathStat.h"

#include <iostream>

int main() { 

	StockPrice sp;
	VSLStreamStatePtr stream = sp.InitGen();

	double *wiener_diff = new double[NSTEPS]; // Random values buffer
	double h = TIME / (double)NSTEPS; // step

	sp.GenerateGauss(0, sqrt(h), NSTEPS, stream, wiener_diff);
	RVCharacteristics rvc(wiener_diff, NSTEPS, h);
	rvc.WriteToCsv(_EPSILON_);
	delete[] wiener_diff;
	sp.FreeGen(stream);

	//StockPrice sp;
	//sp.Execute();
	//std::cout << "Wiener Success" << std::endl;

	//AnSimple as;
	//as.Execute();
	//std::cout << "AnSimple Success" << std::endl;
	//
	//AnExtended ae;
	//ae.Execute();
	//std::cout << "AnExtended Success" << std::endl;

	//NumMethodW nmw;
	//nmw.Execute(nmw.step_array[0], "_Euler_Marayama.csv");
	//std::cout << "Euler Method Success" << std::endl;
	//nmw.Execute(nmw.step_array[1], "_Milstein.csv");
	//std::cout << "Milstein Method Success" << std::endl;
	//nmw.Execute(nmw.step_array[2], "_RK1.csv");
	//std::cout << "RK1 Method Success" << std::endl;

	//NumMethodWZ nmwz;
	//nmwz.Execute(nmwz.step_array[0], "_Burrage_Platen.csv");
	//std::cout << "Burrage Platen Success" << std::endl;
	//nmwz.Execute(nmwz.step_array[1], "_Taylor2.csv");
	//std::cout << "Taylor2 Success" << std::endl;

	system("pause");
	return 0; 
}