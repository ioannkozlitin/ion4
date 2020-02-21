#ifndef MIXDATA_H
#define MIXDATA_H

#include <vector>
#include "../saha/src/elements.h"

class MixData
{
public:
    MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, double T, double V);
    std::vector<TElement> elements;
    std::vector<double> x;
    std::vector<std::vector<double>> xx;
    double T;
    double V;
};

#endif
