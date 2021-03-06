#include "sahamixsolver.h"
#include "sahaleft.h"
#include "findroot.h"

#include <cmath>
#include <limits>

SahaMixResult SahaMixSolver::operator()(const MixData &data)
{
    SahaLeft sahaLeft;
    FindRoot findroot;

    double vFull = data.GetFullV();
    double xe;
    double vFree;

    xe = findroot(-746, log(data.maxZ), [&](double x)
    {
        double vFree = vFreeByXe(data,vFull,x);

        if (vFree <= 0)
        {
            return std::numeric_limits<double>::max();
        }

        return sahaLeft(data, vFree, x);
    }, 1e-7, data.T, vFull);

    vFree = vFreeByXe(data,vFull,xe);

    double vError = fabs((sahaLeft.GetVIon(data, vFree, xe) + vFree)/vFull - 1);

    SahaMixResult rezult1 = {xe, vFree};

    if(vError > 1e-4)
    {
        xe = findroot(-746, log(data.maxZ), xe, [&](double x)
        {
            return xeByVfreeByXe(data, x, vFree) - x;
        }, 1e-7, data.T, vFull);

        xe = xeByVfreeByXe(data, xe, vFree);

        //double vError2 = fabs((sahaLeft.GetVIon(data, vFree, xe) + vFree)/vFull - 1);
        //if(vError2 < vError) return {xe, vFree};
        return {xe, vFree};
    }

    return rezult1;
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

double SahaMixSolver::xeByVfreeByXe(const MixData &data, double xe, double &vFree)
{
    SahaLeft sahaLeft;
    FindRoot findroot;

    double vFull = data.GetFullV();
    vFree = findroot(log(vFull) - 30, log(vFull), [&](double vfree)
    {
        return (sahaLeft.GetVIon(data, vfree, xe) + vfree) / vFull - 1;
    }
    , 1e-12, data.T, vFull);

    double result = sahaLeft(data, vFree, xe);
    return result + xe;
}
