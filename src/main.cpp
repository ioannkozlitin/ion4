#include "mixdata/mixdata.h"
#include "raizermixsolver/raizermixsolver.h"
#include <fstream>
#include <iostream>
#include <cmath>

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

void calculatorMixRaizer(const std::vector<unsigned int> &Z, const std::vector<double> &x, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep, std::string filename)
{
    RaizerMixSolver mixSolver;

    std::vector<double> _lgT;
    std::vector<double> lgVa;
    std::vector<double> _lgRho;

    std::vector<std::vector<double>> ionizationTable;

    for (double lgT = lgTMax; lgT > lgTMin - lgTStep / 2.0; lgT -= lgTStep)
    {
        _lgT.push_back(lgT);
        std::cout << "[" << lgT << "]" << std::flush;
        std::vector<double> ionizationLine;

        bool fillFlag = _lgRho.empty();

        for (double lgRho = lgRhoMax; lgRho > lgRhoMin - lgRhoStep / 2.0; lgRho -= lgRhoStep)
        {
            MixData md(Z, x, pow(10, lgT), pow(10, lgRho));
            ionizationLine.push_back(mixSolver(md));

            if(fillFlag)
            {
                _lgRho.push_back(lgRho);
                lgVa.push_back(log10(md.V));
            }
        }

        ionizationTable.push_back(ionizationLine);
    }

    std::fstream f(filename.c_str(), std::fstream::out);
    f << std::scientific;

    f << "Z=[";for(auto &z : Z) f << z << " ";f << "];" << std::endl;
    f << "x=[";for(auto &_x : x) f << _x << " ";f << "];" << std::endl;

    outputArray(f, "lgT", _lgT);
    outputArray(f, "lgV", lgVa);
    outputArray(f, "lgRho", _lgRho);
    outputTable(f, "xe", ionizationTable);
}

int main()
{	
	try
	{
        calculatorMixRaizer({7, 8, 18}, {0.78, 0.21, 0.01}, -6, 6, 0.1, -2.5, 4.6, 0.1, "result.m");
	}
	catch (std::exception& r)
	{
		printf("\n%s", r.what());
	}
	
    return 0;
}
