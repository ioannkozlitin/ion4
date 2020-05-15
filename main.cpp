#include "saha/src/elements.h"
#include "saha/src/sahasolver.h"
#include "saha/src/atom_ed.h"
#include "saha/saha.h"
#include <math.h>

#include "mix/mixdata.h"
#include "mix/sahamixsolver.h"

#include "saha/src/atom_ed.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void outputArray(std::ostream& os, const std::string dataName, const std::vector<double> &data)
{
   os << dataName << " = [";
   for(std::vector<double>::const_iterator it = data.begin(); it!= data.end(); ++it) {os << *it << " ";};
   os << "];" << std::endl;
}

void outputTable(std::ostream& os, std::string tableName, const std::vector<std::vector<double>> &table)
{
    os << tableName << " = [" << std::endl;
    for (size_t i = 0; i < table.size(); ++i)
    {
        const std::vector<double>& line = table[i];
        for (size_t j = 0; j < line.size(); ++j)
        {
            os << (line[j]) << " ";
        }
        os << std::endl;
    }
    os << "];" << std::endl;
}

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
    const TElement elem(Z, rCoeff, true);

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
    std::vector<std::vector<std::vector<std::vector<double>>>> xxTable;

    xxTable.resize(Z.size());
    for(int i = 0; i < Z.size(); i++) xxTable[i].resize(Z[i] + 1);

    for (double lgT = lgTMax; lgT > lgTMin - lgTStep / 2.0; lgT -= lgTStep)
    {
        lgTPhys.push_back(lgT);
        std::cout << "[" << lgT << "]" << std::flush;
        std::vector<double> ionizationLine;
        std::vector<std::vector<std::vector<double>>> xxLines;

        xxLines.resize(Z.size());
        for(int i = 0; i < Z.size(); i++) xxLines[i].resize(Z[i] + 1);

        bool fillFlag = _lgRho.empty();

        for (double lgRho = lgRhoMax; lgRho > lgRhoMin - lgRhoStep / 2.0; lgRho -= lgRhoStep)
        {
            MixData md(Z, x, rCoeff, true, pow(10, lgT), pow(10, lgRho));
            ionizationLine.push_back(mixSolver.GetFullIonizationInfo(md));

            for(int i = 0; i < Z.size(); i++)
            {
                for(int j = 0; j <= Z[i]; j++) xxLines[i][j].push_back(md.xx[i][j]);
            }

            if(fillFlag)
            {
                _lgRho.push_back(lgRho);
                lgVa.push_back(log10(md.GetFullV()));
            }
        }

        for(int i = 0; i < Z.size(); i++)
        {
            for(int j = 0; j <= Z[i]; j++)
            {
                xxTable[i][j].push_back(xxLines[i][j]);
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

    for(int i = 0; i < Z.size(); i++)
    {
        for(int j = 0; j <= Z[i]; j++)
        {
            outputTable(f, "x_"+std::to_string(Z[i])+"_"+std::to_string(j), xxTable[i][j]);
        }
    }

}

int main()
{	
	try
	{
        calculatorRho_eV(29, 0.6, -6, 6, 0.1, -2.51, 4.6, 0.1, "CuOld.m");
        calculatorMix({29}, {1}, 0.6, -6, 6, 0.1, -2.5, 4.6, 0.1, "CuNew.m");

        //calculatorRho_eV(29, 0.6, 4.199, 4.2, 0.1, -2.51, -2.5, 0.1, "CuTest0.m");
        //calculatorMix({29}, {1}, 0.6, 4.2, 4.2, 0.1, -2.5, -2.5, 0.1, "CuTest1.m");

        //calculatorMix({18,36}, {0.5, 0.5}, 0.6, -6, 6, 0.1, -2.5, 4.6, 0.1, "ArKr.m");
        //calculatorMix({7,8}, {0.79, 0.21}, 0.6, -6, 6, 0.1, 0, 2, 0.01, "air.m");

        //double rho = log10(1.29e-3);
        //calculatorMix({7,8}, {0.79, 0.21}, 0, rho, rho + 0.05, 0.1, 0, 2, 0.01, "air.m");
        //calculatorMix({7,8,18}, {0.7811, 0.2095, 0.0094}, 0, rho, rho + 0.05, 0.1, 0, 2, 0.01, "air.m");
	}
	catch (std::exception& r)
	{
		printf("\n%s", r.what());
	}
	
    return 0;
}
