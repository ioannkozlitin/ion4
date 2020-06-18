#ifndef MIXDATA_H
#define MIXDATA_H

#include <vector>
#include "element.h"

class MixData
{
public:
    MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double T, double V, bool TVaeFlag);
    MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double TeV, double Rho);

    std::vector<Element> elements;
    std::vector<double> x;
    double T;
    double V;
    double maxZ;

private:
    void initZX(const std::vector<unsigned int> &Z, const std::vector<double> &x);
};

#endif // MIXDATA_H
