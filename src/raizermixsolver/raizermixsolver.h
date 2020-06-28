#ifndef RAIZERMIXSOLVER_H
#define RAIZERMIXSOLVER_H

#include "../mixdata/mixdata.h"

class RaizerMixSolver
{
public:
    double operator()(const MixData &data, double eps, bool useBrent, bool useLog);
    int FunCallNum();

private:
    double getF(const MixData &data, double xe);
    double getXePart(const Element &element, double T, double fiPart);
    double getGamma(double fiPart, double fi, double T);
    double getEps(double fi0, double fi1, double T);
    double getXePartOutside(double k, double gamma);
    double getXePartInside(double k, double gamma, double eps);

    int _funCallNum = 0;
};

#endif // RAIZERMIXSOLVER_H
