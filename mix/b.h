#ifndef B_H
#define B_H

#include <vector>
#include "mixdata.h"
#include "dfi.h"

class B
{
public:
    const std::vector<std::vector<double>> &operator()(const MixData &data, double vFree, double xe);

private:
    double p0(double t, double T, double vFree, double xe);

    Dfi dfi;
    std::vector<std::vector<double>> _B;
};

#endif // B_H
