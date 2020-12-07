#ifndef IGMDK_DSRS_H
#define IGMDK_DSRS_H
#include <cmath>
#include "../Utils/Vector.h"
#include "../NumericalMethods/NumericalCommon.h"
#include "../RandomNumberGeneration/Random.h"
namespace igmdk{

template<typename FUNCTION> pair<Vector<double>, double>
    DSRS(Vector<double> const& x0, FUNCTION const& f,
    double step = 1, double factor = 0.8, int maxFEvals = 10000000,
    double yPrecision = numeric_limits<double>::epsilon())
{
    pair<Vector<double>, double> xy(x0, f(x0));
    for(double dd = 0; --maxFEvals && step * (dd + step) > yPrecision;)
    {
        Vector<double> direction = x0;//ensure non-zero direction
        for(int j = 0; j < direction.getSize(); ++j) direction[j] =
            GlobalRNG().uniform01() * GlobalRNG().sign();
        direction *= 1/norm(direction);
        double yNew = f(xy.first + direction * step);
        if(isELess(yNew, xy.second, yPrecision))
        {
            dd = (xy.second - yNew)/step;
            xy.first += direction * step;
            xy.second = yNew;
            step *= 2;
        }
        else step *= factor;
    }
    return xy;
}

}//end namespace igmdk
#endif
