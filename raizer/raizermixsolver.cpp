#include "raizermixsolver.h"
#include "../mix/findroot.h"

#include <cmath>
#include <limits>

double RaizerMixSolver::operator()(const MixData &data)
{
    FindRoot findroot;

    unsigned int maxZ = 0.0;
    for (int j = 0; j < data.elements.size(); j++)
    {
        if (maxZ <= data.elements[j].Z)
        {
            maxZ = data.elements[j].Z;
        }
    }

    double logA = -746;
    double logB = log(maxZ); // log(data.maxZ);

    auto F = [&](double xe)
    {
        return getF(data, xe);
    };

    double eps = 1e-6;
    double T = data.T;
    double V = data.GetFullV();

    double xe = findroot(logA, logB, F, eps, T, V);

    return xe;
}

double RaizerMixSolver::getF(const MixData &data, double xe)
{
    double A = 2 * data.GetFullV() * pow((data.T / (2 * M_PI)), 3.0 / 2.0);
    double fiPart = data.T * log(A / xe);

    double sum = 0.0;
    double xePart;

    for (int i = 0; i < data.x.size(); i++)
    {
        xePart = getXePart(data.elements[i], data.T, fiPart);
        sum += data.x[i] * xePart;
        if (xePart < 0.0)
        {
            printf("\n\nERROR!\n\n"); fflush(0);
        }
    }

    return xe - sum;
}

double RaizerMixSolver::getXePart(const TElement &element, double T, double fiPart)
{
    double gamma, eps, xePart = -1;

    if (fiPart < element.fi[0])
    {
        gamma = exp((fiPart - element.fi[0]) / T);
        xePart = gamma / (gamma + 1.0);
        return xePart;
    }
    else if (fiPart >= element.fi[element.Z - 1])
    {
        gamma = exp((fiPart - element.fi[element.Z - 1]) / T);
        xePart = element.Z - 1.0 / (gamma + 1.0);
        return xePart;
    }
    else
    {
        for (int i = 0; i < element.Z - 1; i++)
        {
            if (element.fi[i] <= fiPart && fiPart < element.fi[i + 1])
            {
                eps = 1.0 / (exp((element.fi[i + 1] - element.fi[i]) / (2.0 * T)) - 1.0);

                if (fiPart < (element.fi[i] + element.fi[i + 1]) / 2.0)
                {
                    gamma = exp((fiPart - element.fi[i]) / T);
                    xePart = (i + 1.0) - (1.0 - (gamma - 1.0) * eps) / (gamma + 1.0);
                }
                else
                {
                    gamma = exp((fiPart - element.fi[i + 1]) / T);
                    xePart = (i + 2.0) - (1.0 - (gamma - 1.0) * eps) / (gamma + 1.0);
                }

                return xePart;
            }
        }
    }
}
