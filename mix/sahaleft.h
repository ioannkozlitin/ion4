#ifndef SAHALEFT_H
#define SAHALEFT_H

#include "mixdata.h"
#include "b.h"

class SahaLeft
{
public:
    double operator()(const MixData &data, double vFree, double xe);
    double GetVIon(const MixData &data, double vFree, double xe);

private:
    B b;
};

#endif // SAHALEFT_H
