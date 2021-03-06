#include "findroot.h"
#include <cmath>
#include <limits>

double FindRoot::operator ()(double logA, double logB, const std::function<double (double)> &F, double eps, double T, double V)
{
    const std::function<double (double)> &Fexp = [&](double x)
    {
        return F(exp(x));
    };

    return core(logA, logB, Fexp(logA), Fexp(logB), Fexp, eps, T, V);
}

double FindRoot::operator ()(double logA, double logB, double x0, const std::function<double (double)> &F, double eps, double T, double V)
{
    const std::function<double (double)> &Fexp = [&](double x)
    {
        return F(exp(x));
    };

    double xprev = x0, fprev = 0;

    if(x0 < 0)
    {
        return core(logA, logB, Fexp(logA), Fexp(logB), Fexp, eps, T, V);
    }

    //return core(logA, logB, F(exp(logA)), F(exp(logB)), log(x0), F, eps, T, V);

    double fx0 = F(x0);

    if(fabs(fx0) < eps) return x0;

    if(!std::isfinite(fx0))
    {
        return core(logA, logB, Fexp(logA), Fexp(logB), Fexp, eps, T, V);
    }

    double x1 = fx0 + x0;
    double fx1 = F(x1);

    if(!std::isfinite(fx1))
    {
        return core(logA, logB, Fexp(logA), Fexp(logB), Fexp, eps, T, V);
    }

    int watchDog = 0;
    while(fabs(fx1) > eps)
    {
        if((((fx0 <= 0) && (fx1 >= 0)) || ((fx0 >= 0) && (fx1 <= 0))) && (x1 > 0))
        {
            if(x1 > x0) return core(log(x0), log(x1), fx0, fx1, Fexp, eps, T, V);
            else return core(log(x1), log(x0), fx1, fx0, Fexp, eps, T, V);
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
                double fa = Fexp(logA);
                double fb = Fexp(logB);

                if(((fa <= 0) && (fx0 >=0)) || ((fa >= 0) && (fx0 <=0)))
                {
                    return core(logA, log(x0), fa, fx0, Fexp, eps, T, V);
                }
                else
                {
                    return core(log(x0), logB, fx0, fb, Fexp, eps, T, V);
                }
            }
        }
    }

    if((fabs(fx1) < eps) && (std::isfinite(fx1))) return x1;
    else
    {
        return core(logA, logB, Fexp(logA), Fexp(logB), Fexp, eps, T, V);
    }
}

double FindRoot::core(double a, double b, double fa, double fb, const std::function<double (double)> &F, double eps, double T, double V)
{
    return core(a, b, fa, fb, 0.5*(a+b), F, eps, T, V);
}

double FindRoot::core(double a, double b, double fa, double fb, double c, const std::function<double (double)> &F, double eps, double T, double V)
{
    return exp(basecore(a,b,fa,fb,c,F,eps,T,V));
}

double FindRoot::basecore(double a, double b, double fa, double fb, double c, const std::function<double (double)> &F, double eps, double T, double V)
{
    double cnew;
    double fc;
    double root = 0;

    int cc = 0, cc2 = 0;

    double Froot;

    cnew = c;
    if(((fa >= 0) && (fb >= 0)) || ((fa <= 0) && (fb <= 0)))
    {
        if (fabs(fa) < fabs(fb))
        {
            root = a;Froot = fa;
        }
        else
        {
            root = b;Froot = fb;
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

            fc = F(c);

            if(cc2 < 30) cnew = invxroot(a,fa,b,fb,c,fc);

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
        while ((b - a > eps) && (cc < 90) && (fabs(fa) > eps) && (fabs(fb) > eps));
        //while ((b - a > eps) && (cc < 90) && (fabs(fc) > eps));
        root = c;Froot = fc;
    }

    if(fabs(Froot) > fabs(fa)) {root = a;Froot = fa;};
    if(fabs(Froot) > fabs(fb)) {root = b;Froot = fb;};

    return root;
}

double FindRoot::fastcore(double a, double b, double fa, double fb, double c, const std::function<double (double)> &F, double eps, double T, double V)
{
    double fc;
    double root = 0;

    int cc = 0, cc2 = 0;

    double x0, x1, x2, x3, fx0, fx1, fx2, fx3;
    double Froot;

    if((c < a) || (c > b)) c = 0.5 * (a + b);

    x0 = a;fx0 = fa;
    x1 = b;fx1 = fb;

    if(((fa >= 0) && (fb >= 0)) || ((fa <= 0) && (fb <= 0)))
    {
        if (fabs(fa) < fabs(fb))
        {
            root = a;Froot = fa;
        }
        else
        {
            root = b;Froot = fb;
        }
    }
    else
    {
        fc = F(c);
        x2 = c;fx2 = fc;
        x0 = a;fx0 = fa;
        x1 = b;fx1 = fb;

        do
        {
            if(cc > 0)
            {
                x3 = invxroot(x0,fx0,x1,fx1,x2,fx2);

                if((x3 > a) && (x3 < b) && (cc2 < 30))
                {
                    fx3 = F(x3);
                    c = x3;fc = fx3;

                    x0 = x1;fx0 = fx1;
                    x1 = x2;fx1 = fx2;
                    x2 = x3;fx2 = fx3;

                    cc2++;
                }
                else
                {
                    c = 0.5 * (a + b);
                    fc = F(c);
                    x2 = c;fx2 = fc;
                    x0 = a;fx0 = fa;
                    x1 = b;fx1 = fb;
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
        while ((b - a > eps) && (cc < 90) && (fabs(fa) > eps) && (fabs(fb) > eps));
        root = c;
        Froot = fc;
    }

    if(fabs(Froot) > fabs(fa)) {root = a;Froot = fa;};
    if(fabs(Froot) > fabs(fb)) {root = b;Froot = fb;};

    return root;
}

double FindRoot::invxroot(double x1, double y1, double x2, double y2,double x3, double y3)
{
    double a=((x3-x2)/(y3-y2)-(x2-x1)/(y2-y1))/(y3-y1);
    double b=a*(y3-y2)+(x3-x2)/(y3-y2);

    return a*y3*y3-b*y3+x3;
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
