#include "raizermixsolver.h"
#include "findroot.h"
#include <cmath>
#include <limits>

double RaizerMixSolver::operator()(const MixData &data)
{
    FindRoot findroot;

    auto F = [&](double xe)
    {
        return getF(data, xe);
    };

    double xe = findroot(-746, log(data.maxZ), F, 1e-6, data.T, data.V);

    return xe;
}

double RaizerMixSolver::getF(const MixData &data, double xe)
{
    double A, fiPart;

    if (xe > std::numeric_limits<double>::min())
    {
        A = 2.0 * data.V * pow((data.T / (2.0 * M_PI)), 3.0 / 2.0);
        fiPart = data.T * log(A / xe);
    }
    else
    {
        fiPart = std::numeric_limits<double>::max();
    }

    double sum = 0.0;
    double xePart;

    for (int i = 0; i < data.x.size(); i++)
    {
        xePart = getXePart(data.elements[i], data.T, fiPart);
        sum += data.x[i] * xePart;

        if (xePart < 0.0)
        {
            printf("\nError: incorrect xe_part\n"); fflush(0);
        }
    }

    return xe - sum;
}

double RaizerMixSolver::getXePart(const Element &element, double T, double fiPart)
{
    double gamma, eps;
    double xePart = -1.0;

    if (fiPart < element.fi[0])
    {
        gamma = getGamma(fiPart, element.fi[0], T);
        xePart = getXePartOutside(1, gamma);
        return xePart;
    }
    else if (fiPart >= element.fi[element.Z - 1])
    {
        gamma = getGamma(fiPart, element.fi[element.Z - 1], T);
        xePart = getXePartOutside(element.Z, gamma);
        return xePart;
    }
    else
    {
        for (int i = 0; i < element.Z - 1; i++)
        {
            if (element.fi[i] <= fiPart && fiPart < element.fi[i + 1])
            {
                eps = getEps(element.fi[i], element.fi[i + 1], T);

                if (fiPart < (element.fi[i] + element.fi[i + 1]) / 2.0)
                {
                    gamma = getGamma(fiPart, element.fi[i], T);
                    xePart = getXePartInside(i + 1, gamma, eps);
                }
                else
                {
                    gamma = getGamma(fiPart, element.fi[i + 1], T);
                    xePart = getXePartInside(i + 2, gamma, eps);
                }

                return xePart;
            }
        }
    }
}

double RaizerMixSolver::getGamma(double fiPart, double fi, double T)
{
    return exp((fiPart - fi) / T);
}

double RaizerMixSolver::getEps(double fi0, double fi1, double T)
{
    return 1.0 / (exp((fi1 - fi0) / (2.0 * T)) - 1.0);
}

double RaizerMixSolver::getXePartOutside(double k, double gamma)
{
    return k - 1.0 / (gamma + 1.0);
}

double RaizerMixSolver::getXePartInside(double k, double gamma, double eps)
{
    return k - (1.0 - (gamma - 1.0) * eps) / (gamma + 1.0);
}
