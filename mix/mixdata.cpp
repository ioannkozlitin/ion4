#include "mixdata.h"
#include <stdio.h>

MixData::MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, double T, double V)
{
    if(x.size() != Z.size())
    {
        printf("Error: size(Z) != size(x)\n");
    }

    for(auto &zi : Z)
    {
        elements.push_back(TElement(zi, rCoeff, correctV0));
        xx.push_back(std::vector<double>(zi+1));
    }

    this->x = x;
    this->T = T;
    this->V = V;
}
