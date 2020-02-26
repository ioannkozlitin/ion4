#include "dfi.h"
#include "../saha/src/mu_fi.h"

#include <cmath>

double Dfi::operator()(const MixData &data, double vFree, double xe, double p0, int i, int j)
{
    const std::vector<double> &v = data.elements[i].v;
    return p0 * (v[j] - v[j-1]);
}
