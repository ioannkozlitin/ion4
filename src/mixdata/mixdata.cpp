#include "mixdata.h"
#include "atomed.h"
#include <stdio.h>

MixData::MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double T, double V, bool TVaeFlag)
{
    if(TVaeFlag)
    {
        this->T = T;
        this->V = V;
    }
    else
    {
        printf("TVaeFlag must be true!\n");
        return;
    }

    initZX(Z, x);
}

MixData::MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double TeV, double Rho)
{
    double meanA = 0;
    for(int i = 0; i < Z.size(); i++)
    {
        meanA += x[i] * GetA(Z[i]);
    }

    this->V = meanA * eRo / Rho;
    this->T = TeV / eFi;

    initZX(Z, x);
}

void MixData::initZX(const std::vector<unsigned int> &Z, const std::vector<double> &x)
{
    if(x.size() != Z.size())
    {
        printf("Error: size(Z) != size(x)\n");
        return;
    }

    for(auto &zi : Z)
    {
        elements.push_back(Element(zi));
    }

    this->x = x;

    maxZ = 0;
    for(int i = 0; i < x.size(); i++) maxZ += Z[i] * x[i];
}
