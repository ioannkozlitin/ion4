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

double calculatorMixRaizer(const std::vector<unsigned int> &Z, const std::vector<double> &x, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep, double eps, bool useBrent, bool useLog, std::vector<std::vector<double>> &ionizationTable, std::string filename)
{
    RaizerMixSolver mixSolver;

    std::vector<double> _lgT;
    std::vector<double> lgVa;
    std::vector<double> _lgRho;

    ionizationTable.clear();

    size_t calcNum = 0;

    for (double lgT = lgTMax; lgT > lgTMin - lgTStep / 2.0; lgT -= lgTStep)
    {
        _lgT.push_back(lgT);
        // std::cout << "[" << lgT << "]" << std::flush;
        std::vector<double> ionizationLine;

        bool fillFlag = _lgRho.empty();

        for (double lgRho = lgRhoMax; lgRho > lgRhoMin - lgRhoStep / 2.0; lgRho -= lgRhoStep)
        {
            MixData md(Z, x, pow(10, lgT), pow(10, lgRho));
            ionizationLine.push_back(mixSolver(md, eps, useBrent, useLog));

            calcNum++;

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

    return mixSolver.FunCallNum() / double(calcNum);
}

double maxDeviation(const std::vector<std::vector<double>> &t1, const std::vector<std::vector<double>> &t2)
{
    double _maxDeviation = -1.0;

    for (size_t i = 0; i < t1.size(); i++)
    {
        for (size_t j = 0; j < t1[i].size(); j++)
        {
            _maxDeviation = std::max(_maxDeviation, fabs(t1[i][j] - t2[i][j]));
        }
    }

    return _maxDeviation;
}

int main()
{	
	try
	{
        std::vector<double> eps = {1e-6};

        for (size_t i = 0; i < eps.size(); i++)
        {
            std::vector<std::vector<double>> ourWithLogXe, brentWithLogXe, ourXe, brentXe;

            double ourWithLogCnt   = calculatorMixRaizer({7, 8, 18}, {0.78, 0.21, 0.01}, -6, 6, 0.1, -2.5, 4.6, 0.1, eps[i], false, true , ourWithLogXe  , "result.m");
            double brentWithLogCnt = calculatorMixRaizer({7, 8, 18}, {0.78, 0.21, 0.01}, -6, 6, 0.1, -2.5, 4.6, 0.1, eps[i], true , true , brentWithLogXe, "result.m");
            double ourCnt          = calculatorMixRaizer({7, 8, 18}, {0.78, 0.21, 0.01}, -6, 6, 0.1, -2.5, 4.6, 0.1, eps[i], false, false, ourXe         , "result.m");
            double brentCnt        = calculatorMixRaizer({7, 8, 18}, {0.78, 0.21, 0.01}, -6, 6, 0.1, -2.5, 4.6, 0.1, eps[i], true , false, brentXe       , "result.m");

            printf("Eps: %g\n"               , eps[i]                                                     );
            printf("Our Log  : %05.2f\n"     , ourWithLogCnt                                              );
            printf("Brent Log: %05.2f | %g\n", brentWithLogCnt, maxDeviation(ourWithLogXe, brentWithLogXe));
            printf("Our      : %05.2f | %g\n", ourCnt         , maxDeviation(ourWithLogXe, ourXe)         );
            printf("Brent    : %05.2f | %g\n", brentCnt       , maxDeviation(ourWithLogXe, brentXe)       );
        }

	}
	catch (std::exception& r)
	{
		printf("\n%s", r.what());
	}
	
    return 0;
}
