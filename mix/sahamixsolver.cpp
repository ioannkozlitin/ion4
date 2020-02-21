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
