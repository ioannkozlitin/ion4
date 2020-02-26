#include "sahamixsolver.h"
#include "sahaleft.h"
#include "findroot.h"

#include <cmath>

double SahaMixSolver::operator()(const MixData &data)
{
    SahaLeft sahaLeft;
    FindRoot findroot;

    double maxZ = 0;
    for(int i = 0; i < data.x.size(); i++) maxZ += data.elements[i].Z * data.x[i];

    double xe = findroot(-746, log(maxZ), [&](double x) {return sahaLeft(data, data.V, x);}, 1e-7, data.T, data.V);

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
