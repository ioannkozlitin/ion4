#include "raizer.h"
#include "raizermixsolver.h"
#include "../mix/sahamixsolver.h"
#include "../output.h"
#include "../saha/src/atom_ed.h"
#include "../saha/src/sahasolver.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <limits>

void calculatorMixRaizer(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep, std::string filename)
{
    RaizerMixSolver mixSolver;

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
            ionizationLine.push_back(mixSolver(md));

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

void getError(const std::vector<double> &x, const std::vector<double> &xRef, double &rms, double &max)
{
    double xMin = std::numeric_limits<double>::max();
    double xMax = -std::numeric_limits<double>::max();

    double xRefMin = std::numeric_limits<double>::max();
    double xRefMax = -std::numeric_limits<double>::max();

    double rmsError = 0.0;

    double curError;
    double maxError = -std::numeric_limits<double>::max();

    for (size_t i = 0; i < x.size(); i++)
    {
        xMin = std::min(xMin, x[i]);
        xMax = std::max(xMax, x[i]);

        xRefMin = std::min(xRefMin, xRef[i]);
        xRefMax = std::max(xRefMax, xRef[i]);

        rmsError += (x[i] - xRef[i]) * (x[i] - xRef[i]);

        curError = std::abs(x[i] - xRef[i]);
        maxError = std::max(maxError, curError);
    }

    double xRange = xMax - xMin;
    double xRefRange = xRefMax - xRefMin;
    double range = std::max(xRange, xRefRange);

    rms = std::sqrt(rmsError / x.size());
    rms = rms / range * 100.0;

    max = maxError / range * 100.0;
}

void raizerTest(int argc, char *argv[])
{
    bool useOldSaha = true;

    double lgTMin  = strtod(argv[1], nullptr);
    double lgTMax  = strtod(argv[2], nullptr);
    double lgTStep = strtod(argv[3], nullptr);

    double lgV = strtod(argv[4], nullptr);

    std::vector<unsigned int> Z;
    std::vector<double> x;

    for (int i = 5; i < argc; i += 2)
    {
        Z.push_back(strtod(argv[i], nullptr));
        x.push_back(strtod(argv[i + 1], nullptr));
    }

    double rCoeff = 0.6;

    SahaMixSolver sahaMixSolver;
    RaizerMixSolver raizerMixSolver;

    std::vector<double> xeSaha, xeRaizer;

    for (double lgT = lgTMax; lgT > lgTMin - lgTStep / 2.0; lgT -= lgTStep)
    {
        MixData mixData(Z, x, rCoeff, true, pow(10, lgT) / eFi, pow(10, lgV), true);

        if (useOldSaha)
        {
            TElement element(Z.front(), rCoeff);
            SahaSolver sahaSolver(element, 0);
            xeSaha.push_back(sahaSolver.Calculate_lgTeV_lgVae(lgT, lgV).Xe);
        }
        else
        {
            xeSaha.push_back(sahaMixSolver(mixData).xe);
        }

        xeRaizer.push_back(raizerMixSolver(mixData));
    }

    double rms, max;
    getError(xeRaizer, xeSaha, rms, max);
    std::cout << std::scientific << "output = [" << rms << " " << max << "]" << std::endl;
}
