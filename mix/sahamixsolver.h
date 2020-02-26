#ifndef SAHAMIXSOLVER_H
#define SAHAMIXSOLVER_H

#include "mixdata.h"

struct SahaMixResult
{
    double xe;
    double Vfree;
};

class SahaMixSolver
{
public:
    SahaMixResult operator()(const MixData &data);
    double GetFullIonizationInfo(MixData &data);
};

#endif // SAHAMIXSOLVER_H
