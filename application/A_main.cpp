
#include "../modules/StockPriece/include/S31_NumMethodsW.h"	
#include "../modules/StockPriece/include/S32_NumMethodsWZ.h"	
#include "../etc/MathStat/EM20_ChiSquared.h"
#include "../etc/RNG/ER21_Normal.h"
#include "../modules/EuropeanOption/AnSolution/O11_CallOption.h"
#include "../modules/EuropeanOption/AnSolution/O12_CallPutOption.h"
#include "../modules/EuropeanOption/NumSolution/O21_MonteCarlo.h"
#include "../modules/EuropeanOption/NumSolution/O22_QuadFormula.h"

#include <iostream>

int main() { 
	//clock_t start, finish;
	//QuadratureFormula qf;
	//qf.Execute();
	//start = clock();
	//for (int i = 0; i < 10000; ++i)
	//	qf.GetRPrice(-5.15f, 6.0f);
	//finish = clock();
	//std::cout << "Time is " << (double)(finish - start) / CLOCKS_PER_SEC << std::endl;


	//std::cout << "Rectangle: "	<< qf.GetRPrice(-5.15f, 6.0f)	<< " Time: " << qf.t << std::endl;
	//std::cout << "Trapeze: "	<< qf.GetTPrice(-5.15f, 6.0f)	<< " Time: " << qf.t << std::endl;
	//std::cout << "Simpson: "	<< qf.GetSPrice(-5.15f, 6.0f)	<< " Time: " << qf.t << std::endl;
	//std::cout << "3/8 Rule: "	<< qf.Get3_8Price(-5.15f, 6.0f) << " Time: " << qf.t << std::endl;

	//MonteCarlo mc;
	//mc.Execute();
	//	std::cout << "Monte-Carlo: " << mc.GetMCPrice(1) << " Time: " << mc.t << std::endl;
	//_sleep(500);
	//CallOption co;
	//co.WriteToCsv(4);
	//CallPutOption cpo;
	//cpo.WriteToCsv(4);
	//float *wiener_diff = new float[NSTEPS]; // Random values buffer
	//NormalGen rng(0, __STEP__, __SEED__);
	//rng.RandomArray(wiener_diff, NSTEPS);
	//ChiSquared cs(wiener_diff, NSTEPS, __STEP__, 35, 0.1f);
	//cs.Execute();
	//delete[] wiener_diff;

	//MCG59 mcg(__SEED__);
	//wiener_diff = new float[NSTEPS]; // Random values buffer
	//mcg.RandomArray(wiener_diff, NSTEPS);
	//RVCharacteristics rvc(wiener_diff, NSTEPS, __STEP__);
	//rvc.WriteToCsv(_EPSILON_);
	////for (int i = 0; i < NSTEPS; ++i)
	////	std::cout << wiener_diff[i] << std::endl;
	//delete[] wiener_diff;

	//float h = TIME / NSTEPS;

	//NormalGen rng(0, h, __SEED__);
	//float *wiener_diff1 = new float[NSTEPS]; // Random values buffer
	//rng.RandomArray(wiener_diff1, NSTEPS);

	//RVCharacteristics rvc1(wiener_diff1, NSTEPS, h);
	//rvc1.WriteToCsv(_EPSILON_);

	//_sleep(1000);
	//StockPrice sp;
	//VSLStreamStatePtr stream = sp.InitGen();
	//float *wiener_diff2 = new float[NSTEPS]; // Random values buffer
	//sp.GenerateGauss(0, sqrt(h), NSTEPS, stream, wiener_diff2);
	//RVCharacteristics rvc2(wiener_diff2, NSTEPS, h);
	//rvc2.WriteToCsv(_EPSILON_);

	//for (int i = 0; i < NSTEPS; ++i)
	//	std::cout << wiener_diff[i] << "\t" << wiener_diff1[i] << "\t" << wiener_diff2[i] << std::endl;
	//delete[] wiener_diff1;
	//delete[] wiener_diff;
	//delete[] wiener_diff2;

	//std::cout << "Wiener Success" << std::endl;

	//AnSimple as;
	//as.Execute();
	//std::cout << "AnSimple Success" << std::endl;

	//AnExtended ae;
	//ae.Execute();
	//std::cout << "AnExtended Success" << std::endl;

	NumMethodW nmw;
	nmw.Execute(nmw.StepArray[0], "_Euler_Marayama.csv");
	std::cout << "Euler Method Finished" << std::endl << std::endl;
	nmw.Execute(nmw.StepArray[1], "_Milstein.csv");
	std::cout << "Milstein Method Finished" << std::endl << std::endl;
	nmw.Execute(nmw.StepArray[2], "_RK1.csv");
	std::cout << "RK1 Method Finished" << std::endl << std::endl;

	//NumMethodWZ nmwz;
	//nmwz.Execute(nmwz.step_array[0], "_Burrage_Platen.csv");
	//std::cout << "Burrage Platen Success" << std::endl;
	//nmwz.Execute(nmwz.step_array[1], "_Taylor2.csv");
	//std::cout << "Taylor2 Success" << std::endl;

	system("pause");
	return 0; 
}