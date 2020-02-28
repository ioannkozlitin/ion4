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

private:
    double vFreeByXe(const MixData &data, double V, double xe);
    SahaMixResult resultForVfree(const MixData &data, double vFree, double maxZ);
};

#endif // SAHAMIXSOLVER_H
