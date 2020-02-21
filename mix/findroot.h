#ifndef FINDROOT_H
#define FINDROOT_H

#include <functional>

class FindRoot
{
public:
    double operator()(double logA, double logB, const std::function<double(double)> &F, double eps, double T, double V);

private:
    double xroot(double x1, double y1, double x2, double y2,double x3, double y3);
    double chord(double x1, double y1, double x2, double y2);
};

#endif // FINDROOT_H
