#include "saha/src/elements.h"
#include "saha/src/sahasolver.h"
#include "saha/src/atom_ed.h"
#include "saha/saha.h"
#include <math.h>

#include "mix/mixdata.h"
#include "mix/sahamixsolver.h"
#include "raizer/raizer.h"
#include "output.h"

#include "saha/src/atom_ed.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void CrashTest(double rCoeff, double lgVMin, double lgVMax, double lgVStep, double lgTMin, double lgTMax, double lgTStep)
{
	for (int Z = 1; Z <= 103; Z++)
	{
		const TElement elem(Z, rCoeff);
        SahaSolver solver(elem, 0);

		printf("[%d] ",Z);

		for (double lgT = lgTMax; lgT > lgTMin; lgT -= lgTStep)
		{
			for (double lgV = lgVMin; lgV < lgVMax; lgV += lgVStep)
			{
				SahaPoint res = solver.Calculate_lgTeV_lgVae(lgT, lgV);
			}
		}

	}
}

void calculator(unsigned int Z, double rCoeff, double lgVMin, double lgVMax, double lgVStep, double lgTMin, double lgTMax, double lgTStep, std::string filename)
{
    const TElement elem(Z, rCoeff);

    std::vector<double> lgTPhys;
    std::vector<double> lgVa;
    std::vector<double> lgRho;

    std::vector<std::vector<double>> ionizationTable;
    std::vector<std::vector<double>> pTable;
    std::vector<std::vector<double>> eTable;

    for (double lgT = lgTMax; lgT > lgTMin; lgT -= lgTStep) lgTPhys.push_back(lgT);

    double roConst = log10(eRo*elem.A);
    for (double lgV = lgVMin; lgV < lgVMax; lgV += lgVStep)
    {
        lgVa.push_back(lgV);
        lgRho.push_back(roConst - lgV);
    }

    for (double lgT = lgTMax; lgT > lgTMin; lgT -= lgTStep)
    {
        std::cout << "[" << lgT << "]" << std::flush;
        std::vector<double> ionizationLine;
        std::vector<double> pLine;
        std::vector<double> eLine;
        for (double lgV = lgVMin; lgV < lgVMax; lgV += lgVStep)
        {
            TElement elemCopy = elem;
            //elemCopy.softIon(pow(10,lgV));
            SahaSolver solver(elemCopy, 0);
            SahaPoint res = solver.Calculate_lgTeV_lgVae(lgT,lgV);
            ionizationLine.push_back(res.Xe);
            pLine.push_back(res.P);
            eLine.push_back(res.E);
        }
        ionizationTable.push_back(ionizationLine);
        pTable.push_back(pLine);
        eTable.push_back(eLine);
    }

    std::fstream f(filename.c_str(), std::fstream::out);
    f << std::scientific;

    f << "Z=" << Z << ";" << std::endl;
    outputArray(f, "lgT", lgTPhys);
    outputArray(f, "lgV", lgVa);
    outputArray(f, "lgRho", lgRho);
    outputTable(f, "xe_Saha", ionizationTable);
    outputTable(f, "P_Saha", pTable);
    outputTable(f, "E_Saha", eTable);
}

void calculatorRho_eV(unsigned int Z, double rCoeff, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep, std::string filename)
{
    double roConst = log10(eRo*elements::GetA(Z));
    double lgVMin = roConst - lgRhoMax;
    double lgVMax = roConst - lgRhoMin;
    double lgVStep = lgRhoStep;

    calculator(Z, rCoeff, lgVMin, lgVMax, lgVStep, lgTMin, lgTMax, lgTStep, filename);
}

void calculatorMix(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep, std::string filename)
{
    SahaMixSolver mixSolver;

    std::vector<double> lgTPhys;
    std::vector<double> lgVa;
    std::vector<double> _lgRho;

    std::vector<std::vector<double>> ionizationTable;

    for (double lgT = lgTMax; lgT > lgTMin - lgTStep / 2.0; lgT -= lgTStep)
    {
        lgTPhys.push_back(lgT);
        std::cout << "[" << lgT << "]" << std::flush;
        std::vector<double> ionizationLine;

        bool fillFlag = _lgRho.empty();

        for (double lgRho = lgRhoMax; lgRho > lgRhoMin - lgRhoStep / 2.0; lgRho -= lgRhoStep)
        {
            MixData md(Z, x, rCoeff, true, pow(10, lgT), pow(10, lgRho));
            ionizationLine.push_back(mixSolver(md).xe);

            if(fillFlag)
            {
                _lgRho.push_back(lgRho);
                lgVa.push_back(log10(md.GetFullV()));
            }
        }

        ionizationTable.push_back(ionizationLine);
    }

    std::fstream f(filename.c_str(), std::fstream::out);
    f << std::scientific;

    f << "Z=[";for(auto &z : Z) f << z << " ";f << "];" << std::endl;
    f << "x=[";for(auto &_x : x) f << _x << " ";f << "];" << std::endl;

    outputArray(f, "lgT", lgTPhys);
    outputArray(f, "lgV", lgVa);
    outputArray(f, "lgRho", _lgRho);
    outputTable(f, "xe_Saha", ionizationTable);
}

int main(int argc, char *argv[])
{	  
    // raizerTest(argc, argv); return 0;

	try
	{
        /*SahaMixSolver sms;

        printf("xe0 = [");

        double lgT = 1;

        for(double lgRho = -3; lgRho < 5.0005; lgRho += 0.01)
        {
            //MixData md({18, 36}, {0.5, 0.5}, 0.6, true, pow(10, lgT), 1);
            MixData md({18, 36}, {0.5, 0.5}, 0.6, true, pow(10, lgT), pow(10, lgRho));
            double xe = sms.GetFullIonizationInfo(md);
            printf("%g ", xe);
        }
        printf("];\n");*/

        //calculatorRho_eV(29, 0.6, -6, 6, 0.1, -2.501, 4.6, 0.1, "Cu.m");
        //calculatorMix({29}, {1}, 0.6, -6, 6, 0.1, -2.5, 4.6, 0.1, "CuNew.m");
        calculatorMixRaizer({18, 36}, {0.5, 0.5}, 0.6, -6, 6, 0.1, -2.5, 4.6, 0.1, "Raizer.m");
        //calculatorMix({18,36}, {0.5, 0.5}, 0.6, -6, 6, 0.1, -2.5, 4.6, 0.1, "ArKr.m");
	}
	catch (std::exception& r)
	{
		printf("\n%s", r.what());
	}
	
    return 0;
}
