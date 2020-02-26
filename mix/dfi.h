#ifndef DFI_H
#define DFI_H

#include <vector>
#include "mixdata.h"

class Dfi
{
public:
    double operator()(const MixData &data, double vFree, double xe, double p0, int i, int j);
};

#endif // DFI_H
