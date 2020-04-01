#ifndef RAIZERMIXSOLVER_H
#define RAIZERMIXSOLVER_H

#include "../mix/mixdata.h"

class RaizerMixSolver
{
public:
    double operator()(const MixData &data);

private:
    double getF(const MixData &data, double xe);
    double getXePart(const TElement &element, double T, double fiPart);
};

#endif // RAIZERMIXSOLVER_H
