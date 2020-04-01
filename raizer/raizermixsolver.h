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
    double getGamma(double fiPart, double fi, double T);
    double getEps(double fi0, double fi1, double T);
    double getXePartOutside(double k, double gamma);
    double getXePartInside(double k, double gamma, double eps);
};

#endif // RAIZERMIXSOLVER_H
