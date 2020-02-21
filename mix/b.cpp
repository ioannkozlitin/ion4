#include "b.h"
#include "../saha/src/mu_fi.h"

const std::vector<std::vector<double>> &B::operator()(const MixData &data, double vFree, double xe)
{
    double _mu = mu(data.T, vFree, xe);
    _B.resize(data.xx.size());

    for(int i = 0; i < data.xx.size(); i++)
    {
        _B[i].resize(data.xx[i].size());
        _B[i][0] = 0;
        for(int j = 1; j < data.xx[i].size(); j++)
        {
            double A = data.elements[i].logG[j] - data.elements[i].logG[j-1] - (_mu + data.elements[i].fi[j-1] + dfi(data, vFree, xe, i, j-1)) / data.T;
            _B[i][j] = _B[i][j-1] + A;
        }
    }

    return _B;
}
