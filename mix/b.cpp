#include "b.h"
#include "../saha/src/mu_fi.h"
#include <cmath>

const std::vector<std::vector<double>> &B::operator()(const MixData &data, double vFree, double xe)
{
    double _mu = mu(data.t, vFree, xe);
    _B.resize(data.xx.size());

    double _p0 = p0(data.t, vFree, xe);

    for(int i = 0; i < data.xx.size(); i++)
    {
        _B[i].resize(data.xx[i].size());
        _B[i][0] = 0;
        for(int j = 1; j < data.xx[i].size(); j++)
        {
            double A = data.elements[i].logG[j] - data.elements[i].logG[j-1] - (_mu + data.elements[i].fi[j-1] + dfi(data, vFree, xe, _p0, i, j)) / data.T;
            _B[i][j] = _B[i][j-1] + A;
        }
    }

    return _B;
}

double B::p0(double t, double vFree, double xe)
{
    return 2*sqrt(2.0)/(3*M_PI*M_PI) * pow(t,2.5) * I15mu_d_t(t,vFree,xe) + t / vFree;
}
