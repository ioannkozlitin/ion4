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
    double xeByVfreeByXe(const MixData &data, double xe, double &vFree);
};

#endif // SAHAMIXSOLVER_H
