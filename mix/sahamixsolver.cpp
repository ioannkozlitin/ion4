#include "sahamixsolver.h"
#include "sahaleft.h"
#include "findroot.h"

#include <cmath>

double SahaMixSolver::operator()(const MixData &data)
{
    SahaLeft sahaLeft;
    FindRoot findroot;

    unsigned int maxZ = data.elements[0].Z;
    for(int i = 1; i < data.x.size(); i++) if(data.elements[i].Z > maxZ) maxZ = data.elements[i].Z;

    double xe = findroot(-746, log(double(maxZ)), [&](double x) {return sahaLeft(data, data.V, x);}, 1e-7, data.T, data.V);

    return xe;
}

double SahaMixSolver::GetFullIonizationInfo(MixData &data)
{
    double xe = operator ()(data);

    B b;
    const std::vector<std::vector<double>> &bb = b(data, data.V, xe);

    for(int i = 0; i < bb.size(); i++)
    {
        double siBfull = 0;
        for(int j = 0; j < bb[i].size(); j++)
        {
            siBfull += exp(bb[i][j]);
        }

        for(int j = 0; j < bb[i].size(); j++)
        {
            data.xx[i][j] = data.x[i] / siBfull * exp(bb[i][j]);
        }
    }

    return xe;
}
