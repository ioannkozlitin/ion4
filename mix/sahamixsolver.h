#ifndef SAHAMIXSOLVER_H
#define SAHAMIXSOLVER_H

#include "mixdata.h"

class SahaMixSolver
{
public:
    double operator()(const MixData &data);
    double GetFullIonizationInfo(MixData &data);
};

#endif // SAHAMIXSOLVER_H
