#include "../modules/StockPriece/include/S10_Wiener.h"
#include "../modules/StockPriece/include/S20_AnSimple.h"
#include "../modules/StockPriece/include/S21_AnExtended.h"
#include "../modules/StockPriece/include/S31_NumMethodsW.h"	
#include "../modules/StockPriece/include/S32_NumMethodsWZ.h"	
#include "../etc/MathStat/EM10_MathStat.h"
#include "../etc/RNG/ER10_RNG.h"

#include <iostream>

int main() { 

	//MCG59 mcg(__SEED__);
	//float *wiener_diff = new float[NSTEPS]; // Random values buffer
	//float h = TIME / (float)NSTEPS; // step
	//mcg.RandomArray(wiener_diff, NSTEPS);
	//RVCharacteristics rvc(wiener_diff, NSTEPS, h);
	//rvc.WriteToCsv(_EPSILON_);
	//for (int i = 0; i < NSTEPS; ++i)
	//	std::cout << wiener_diff[i] << std::endl;
	//delete[] wiener_diff;



	//StockPrice sp;
	//VSLStreamStatePtr stream = sp.InitGen();
	//float *wiener_diff = new float[NSTEPS]; // Random values buffer
	//float h = TIME / (float)NSTEPS; // step
	//sp.GenerateGauss(0, sqrt(h), NSTEPS, stream, wiener_diff);
	//RVCharacteristics rvc(wiener_diff, NSTEPS, h);
	//rvc.WriteToCsv(_EPSILON_);
	//delete[] wiener_diff;
	//sp.FreeGen(stream);

	//StockPrice sp1;
	//sp1.Execute();
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