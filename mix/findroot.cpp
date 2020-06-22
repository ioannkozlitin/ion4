#include "findroot.h"
#include <cmath>
#include <limits>

double FindRoot::operator ()(double logA, double logB, const std::function<double (double)> &F, double eps, double T, double V)
{
    //return core0(logA, logB, F(exp(logA)), F(exp(logB)), F, eps, T, V);
    return core(logA, logB, F(exp(logA)), F(exp(logB)), F, eps, T, V);
}

double FindRoot::operator ()(double logA, double logB, double x0, const std::function<double (double)> &F, double eps, double T, double V)
{
    double xprev = x0, fprev = 0;

    if(x0 < 0)
    {
        return core(logA, logB, F(exp(logA)), F(exp(logB)), F, eps, T, V);
    }

    double fx0 = F(x0);

    if(fabs(fx0) < eps) return x0;

    if(!std::isfinite(fx0))
    {
        return core(logA, logB, F(exp(logA)), F(exp(logB)), F, eps, T, V);
    }

    double x1 = fx0 + x0;
    double fx1 = F(x1);

    if(!std::isfinite(fx1))
    {
        return core(logA, logB, F(exp(logA)), F(exp(logB)), F, eps, T, V);
    }

    int watchDog = 0;
    while(fabs(fx1) > eps)
    {
        if((((fx0 <= 0) && (fx1 >= 0)) || ((fx0 >= 0) && (fx1 <= 0))) && (x1 > 0))
        {
            if(x1 > x0) return core(log(x0), log(x1), fx0, fx1, F, eps, T, V);
            else return core(log(x1), log(x0), fx1, fx0, F, eps, T, V);
        }
        else
        {
            if((fabs(fx1) < fabs(fx0)) && (watchDog < 100))
            {
                double fx2, x2;
                if(watchDog < 1)
                {
                    x2 = fx1 + x1;
                    fx2 = F(x2);
                }
                else
                {
                    x2 = xroot(xprev, fprev, x0, fx0, x1, fx1);

                    if(std::isfinite(x2)) fx2 = F(x2);
                    else
                    {
                        x2 = fx1 + x1;
                        fx2 = F(x2);
                    }
                }

                xprev = x0;
                fprev = fx0;
                x0 = x1;
                fx0 = fx1;
                x1 = x2;
                fx1 = fx2;

                if(!std::isfinite(fx1)) break;

                watchDog++;
            }
            else
            {
                double fa = F(exp(logA));
                double fb = F(exp(logB));

                if(((fa <= 0) && (fx0 >=0)) || ((fa >= 0) && (fx0 <=0)))
                {
                    return core(logA, log(x0), fa, fx0, F, eps, T, V);
                }
                else
                {
                    return core(log(x0), logB, fx0, fb, F, eps, T, V);
                }
            }
        }
    }

    if((fabs(fx1) < eps) && (std::isfinite(fx1))) return x1;
    else
    {
        return core(logA, logB, F(exp(logA)), F(exp(logB)), F, eps, T, V);
    }
}

double FindRoot::core0(double a, double b, double fa, double fb, const std::function<double (double)> &F, double eps, double T, double V)
{
    double c, cnew;
    double fc;
    double root = 0;

    int cc = 0, cc2 = 0;

    double Froot;

    cnew = b + 1;
    if(((fa >= 0) && (fb >= 0)) || ((fa <= 0) && (fb <= 0)))
    {
        if (fabs(fa) < fabs(fb))
        {
            root = exp(a);Froot = fa;
        }
        else
        {
            root = exp(b);Froot = fb;
        }
    }
    else
    {
        do
        {
            if((cnew > a) && (cnew < b))
            {
                c = cnew;
                cc2++;
            }
            else
            {
                c = 0.5*(a + b);
            }

            fc = F(exp(c));

            if(cc2 < 16) cnew = xroot(a,fa,b,fb,c,fc);

            if (((fa <= 0) && (fc >= 0)) || ((fa >= 0) && (fc <= 0)))
            {
                b = c; fb = fc;
            }
            else if (((fb <= 0) && (fc >= 0)) || ((fb >= 0) && (fc <= 0)))
            {
                a = c; fa = fc;
            }
            else
            {
                printf("find root error: ln(a) = %g ln(b) = %g ln(c) = %g fa = %g fb = %g fc = %g T = %g V = %g\n", a, b, c, fa, fb, fc, T, V);
                return 0;
            }

            cc++;
        }
        while ((b - a > eps) && (cc < 60) && (fabs(fa) > eps) && (fabs(fb) > eps));
        root = exp(c);Froot = fc;
    }

    if(fabs(Froot) > fabs(fa)) {root = exp(a);Froot = fa;};
    if(fabs(Froot) > fabs(fb)) {root = exp(b);Froot = fb;};

    double root2 = exp(chord(a,fa,b,fb));
    if(fabs(F(root2)) < fabs(Froot))
    {
        return root2;
    }
    else return root;
}

double FindRoot::core(double a, double b, double fa, double fb, const std::function<double (double)> &F, double eps, double T, double V)
{
    double c, cnew;
    double fc;
    double root = 0;

    int cc = 0, cc2 = 0;
    int size = 2;

    double x0, x1, x2, x3, fx0, fx1, fx2, fx3;
    double Froot;

    x0 = a;fx0 = fa;
    x1 = b;fx1 = fb;

    cnew = b + 1;
    if(((fa >= 0) && (fb >= 0)) || ((fa <= 0) && (fb <= 0)))
    {
        if (fabs(fa) < fabs(fb))
        {
            root = exp(a);Froot = fa;
        }
        else
        {
            root = exp(b);Froot = fb;
        }
    }
    else
    {
        do
        {
            if((size < 3))
            {
                c = 0.5 * (a + b);
                fc = F(exp(c));

                size = 3;
                x2 = c;fx2 = fc;
                x0 = a;fx0 = fa;
                x1 = b;fx1 = fb;
            }
            else
            {
                x3 = xroot(x0,fx0,x1,fx1,x2,fx2);

                if((x3 > a) && (x3 < b))
                {
                    fx3 = F(exp(x3));
                    c = x3;fc = fx3;

                    x0 = x1;fx0 = fx1;
                    x1 = x2;fx1 = fx2;
                    x2 = x3;fx2 = fx3;

                    cc2++;
                }
                else
                {
                    //printf("<%d %d %g(%g) %g(%g) %g(%g)>", cc, cc2, a, fa, b, fb, c, fc);
                    size = 2;continue;
                }
            }

            if (((fa <= 0) && (fc >= 0)) || ((fa >= 0) && (fc <= 0)))
            {
                b = c; fb = fc;
            }
            else if (((fb <= 0) && (fc >= 0)) || ((fb >= 0) && (fc <= 0)))
            {
                a = c; fa = fc;
            }
            else
            {
                printf("find root error: ln(a) = %g ln(b) = %g ln(c) = %g fa = %g fb = %g fc = %g T = %g V = %g\n", a, b, c, fa, fb, fc, T, V);
                return 0;
            }

            cc++;
        }
        while ((b - a > eps) && (cc < 60) && (fabs(fa) > eps) && (fabs(fb) > eps));
        root = exp(c);
        Froot = fc;
    }

    if(fabs(Froot) > fabs(fa)) {root = exp(a);Froot = fa;};
    if(fabs(Froot) > fabs(fb)) {root = exp(b);Froot = fb;};

    /*
    double root2 = exp(chord(a,fa,b,fb));

    double F0 = F(root), F2 = F(root2);
    if(fabs(F2) < fabs(F0))
    {
        //printf("[%g %g]", root, root2);
        return root2;
    }
    else return root;
    */

    return root;
}

double FindRoot::xroot(double x1, double y1, double x2, double y2,double x3, double y3)
{
    double a=((y3-y2)/(x3-x2)-(y2-y1)/(x2-x1))/(x3-x1);
    double b=a*(x3-x2)+(y3-y2)/(x3-x2);
    double c=y3;
    double d = b*b-4*a*c;
    double r1 = (-b + sqrt(d))/(2*a);
    double r2 = (-b - sqrt(d))/(2*a);

    if(fabs(r1)<fabs(r2)) return r1+x3;
    else return r2+x3;
}

double FindRoot::chord(double x1, double y1, double x2, double y2)
{
    return x1 - (x2-x1) / (y2-y1) * y1;
}
