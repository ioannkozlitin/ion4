#include "sahamixsolver.h"
#include "sahaleft.h"
#include "findroot.h"

#include <cmath>
#include <limits>

SahaMixResult SahaMixSolver::operator()(const MixData &data)
{
    //printf("[%g %g]", vFreeByXe(data, data.GetFullV(), 27), data.GetFullV());

    SahaLeft sahaLeft;
    FindRoot findroot;

    double vFull = data.GetFullV();
    double xe = findroot(-746, log(data.maxZ), [&](double x)
    {
        double vFree = vFreeByXe(data,vFull,x);

        if (vFree <= 0)
        {
            return std::numeric_limits<double>::max();
        }

        return sahaLeft(data, vFree, x);
    }, 1e-7, data.T, vFull);

    double vFree = vFreeByXe(data,vFull,xe);
    //printf("(%g)", (sahaLeft.GetVIon(data, vFree, xe) + vFree - vFull)/vFull);

    /*
    double vFree = findroot(log(vFull) - 30, log(vFull), [&](double vfree) {return (vfree - resultForVfree(data,vfree,maxZ).Vfree) / vFull;},1e-12, data.T, vFull);
    SahaMixResult result = resultForVfree(data, vFree, maxZ);

    printf("(%g)", (vFree-result.Vfree)/vFull);
    */

    //printf("(%g)", (vFull - vIon) / vFull);

    return {xe, vFree};
}

double SahaMixSolver::GetFullIonizationInfo(MixData &data)
{
    SahaMixResult result = operator ()(data);

    B b;
    const std::vector<std::vector<double>> &bb = b(data, result.Vfree, result.xe);

    for(int i = 0; i < bb.size(); i++)
    {
        double biMax = bb[i][0];
        for(int j = 1; j < bb[i].size(); j++) if(biMax < bb[i][j]) biMax = bb[i][j];

        double siBfull = 0;
        for(int j = 0; j < bb[i].size(); j++)
        {
            siBfull += exp(bb[i][j]-biMax);
        }

        for(int j = 0; j < bb[i].size(); j++)
        {
            data.xx[i][j] = data.x[i] / siBfull * exp(bb[i][j]-biMax);
        }
    }

    return result.xe;
}

double SahaMixSolver::vFreeByXe(const MixData &data, double V, double xe)
{
    double vFree = 0;

    for(int i = 0; i < data.x.size(); i++)
    {
        double xei = xe * data.elements[i].Z / data.maxZ;

        unsigned int intXe = floor(xei);
        double fracXe = xei - intXe;

        if(intXe < data.elements[i].Z)
        {
            vFree += (V - (data.elements[i].v[intXe] * (1-fracXe) + data.elements[i].v[intXe+1] * fracXe)) * data.x[i];
        }
        else
        {
            vFree += (V - data.elements[i].v[data.elements[i].Z]) * data.x[i];
        }
    }

    return vFree;
}

SahaMixResult SahaMixSolver::resultForVfree(const MixData &data, double vFree, double maxZ)
{
    SahaLeft sahaLeft;
    FindRoot findroot;

    double xe = findroot(-746, log(maxZ), [&](double x) {return sahaLeft(data, vFree, x);}, 1e-7, data.T, vFree);
    double vIon = sahaLeft.GetVIon(data, vFree, xe);
    return {xe, data.GetFullV()-vIon};
}
