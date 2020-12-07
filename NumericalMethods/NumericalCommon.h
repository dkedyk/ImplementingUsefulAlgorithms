#ifndef IGMDK_NUMERICAL_COMMON_H
#define IGMDK_NUMERICAL_COMMON_H

#include <cmath>
namespace igmdk{

double defaultPrecEps = sqrt(numeric_limits<double>::epsilon());
double highPrecEps = 100 * numeric_limits<double>::epsilon();

bool isELess(double a, double b,
    double eRelAbs = numeric_limits<double>::epsilon())
    {return a < b && b - a >= eRelAbs * max(1.0, max(abs(a), abs(b)));}
bool isEEqual(double a, double b,
    double eRelAbs = numeric_limits<double>::epsilon())
    {return !isELess(a, b, eRelAbs) && !isELess(b, a, eRelAbs);}

template<typename X> double normInf(Vector<X> const& x)
{//works for complex vector too
    double xInf = 0;
    for(int i = 0; i < x.getSize(); ++i)
    {
        double ax = abs(x[i]);
        if(isnan(ax)) return ax;//check for NaN before max
        xInf = max(xInf, ax);
    }
    return xInf;
}

}//end namespace
#endif
