#include "mixdata.h"
#include <stdio.h>
#include "../saha/src/atom_ed.h"

MixData::MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, double T, double V, bool TVaeFlag)
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

    initZX(Z, x, rCoeff, correctV0);
}

MixData::MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, double TeV, double Rho)
{
    double meanA = 0;
    for(int i = 0; i < Z.size(); i++)
    {
        meanA += x[i] * elements::GetA(Z[i]);
    }

    this->V = meanA * eRo / Rho;
    this->T = TeV / eFi;

    initZX(Z, x, rCoeff, correctV0);
}

void MixData::initZX(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0)
{
    if(x.size() != Z.size())
    {
        printf("Error: size(Z) != size(x)\n");
        return;
    }

    for(auto &zi : Z)
    {
        elements.push_back(TElement(zi, rCoeff, correctV0));
        xx.push_back(std::vector<double>(zi+1));
    }

    this->x = x;

    maxZ = 0;
    for(int i = 0; i < x.size(); i++) maxZ += Z[i] * x[i];

    //Soft Ions
    //for(TElement &elem : elements) elem.softIon(V);
}
